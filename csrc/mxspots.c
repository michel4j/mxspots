#define MXSPOTS_EXPORTS
#include "mxspots.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define TILE_SIZE 32
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int mxspots_get_version(void) {
    return 100; /* Version 1.0.0 encoded as 100 */
}

int mxspots_ping(const MxSpotsParams *params) {
    if (params == NULL) {
        return -1;
    }
    if (params->snr_threshold <= 0.0f) {
        return -2;
    }
    return 0;
}

static int compare_spots_desc(const void *a, const void *b) {
    const MxSpot *sa = (const MxSpot *)a;
    const MxSpot *sb = (const MxSpot *)b;
    if (sb->intensity > sa->intensity) return 1;
    if (sb->intensity < sa->intensity) return -1;
    return 0;
}

int mxspots_find_spots(
    const float *data,
    int nx,
    int ny,
    const MxSpotsParams *params,
    MxSpot *out_spots,
    int max_spots
) {
    if (data == NULL || nx <= 0 || ny <= 0 || params == NULL) {
        return 0;
    }

    int tx_count = (nx + TILE_SIZE - 1) / TILE_SIZE;
    int ty_count = (ny + TILE_SIZE - 1) / TILE_SIZE;
    int total_tiles = tx_count * ty_count;

    float *tile_mean = (float *)malloc(total_tiles * sizeof(float));
    float *tile_std = (float *)malloc(total_tiles * sizeof(float));
    if (tile_mean == NULL || tile_std == NULL) {
        free(tile_mean);
        free(tile_std);
        return 0;
    }

    /* Pass 1 & 2: Multithreaded tile background mean and std calculation */
    #pragma omp parallel for schedule(static)
    for (int ty = 0; ty < ty_count; ++ty) {
        int y_start = ty * TILE_SIZE;
        int y_end = (y_start + TILE_SIZE < ny) ? (y_start + TILE_SIZE) : ny;

        for (int tx = 0; tx < tx_count; ++tx) {
            int x_start = tx * TILE_SIZE;
            int x_end = (x_start + TILE_SIZE < nx) ? (x_start + TILE_SIZE) : nx;
            int t_idx = ty * tx_count + tx;

            double sum = 0.0;
            double sum_sq = 0.0;
            int count = 0;

            for (int y = y_start; y < y_end; ++y) {
                const float *row = &data[y * nx];
                for (int x = x_start; x < x_end; ++x) {
                    float v = row[x];
                    if (v >= 0.0f && !isnan(v) && !isinf(v)) {
                        sum += v;
                        sum_sq += (double)v * v;
                        count++;
                    }
                }
            }

            if (count > 0) {
                double mean1 = sum / count;
                double var1 = (sum_sq / count) - (mean1 * mean1);
                double std1 = (var1 > 0.0) ? sqrt(var1) : 0.0;

                /* Pass 2: Outlier rejection for robust background */
                double sum2 = 0.0;
                double sum_sq2 = 0.0;
                int count2 = 0;
                double cutoff = mean1 + 3.0 * std1;

                for (int y = y_start; y < y_end; ++y) {
                    const float *row = &data[y * nx];
                    for (int x = x_start; x < x_end; ++x) {
                        float v = row[x];
                        if (v >= 0.0f && v <= cutoff && !isnan(v) && !isinf(v)) {
                            sum2 += v;
                            sum_sq2 += (double)v * v;
                            count2++;
                        }
                    }
                }

                if (count2 > 0) {
                    double mean2 = sum2 / count2;
                    double var2 = (sum_sq2 / count2) - (mean2 * mean2);
                    double std2 = (var2 > 0.0) ? sqrt(var2) : 1.0;
                    tile_mean[t_idx] = (float)mean2;
                    tile_std[t_idx] = (std2 < 1.0) ? 1.0f : (float)std2;
                } else {
                    tile_mean[t_idx] = (float)mean1;
                    tile_std[t_idx] = (std1 < 1.0) ? 1.0f : (float)std1;
                }
            } else {
                tile_mean[t_idx] = 0.0f;
                tile_std[t_idx] = 1.0f;
            }
        }
    }

    /* Compute resolution radial limits */
    float distance = (params->distance > 0.0f) ? params->distance : 100.0f;
    float wavelength = params->wavelength;
    float r_min_sq = 0.0f;
    float r_max_sq = 1e30f;

    if (wavelength > 0.0f && distance > 0.0f) {
        if (params->d_max > 0.0f) {
            float s_low = wavelength / (2.0f * params->d_max);
            if (s_low > 0.0f && s_low < 1.0f) {
                float theta_low = asinf(s_low);
                float r_low = distance * tanf(2.0f * theta_low);
                r_min_sq = r_low * r_low;
            }
        }
        if (params->d_min > 0.0f) {
            float s_high = wavelength / (2.0f * params->d_min);
            if (s_high > 0.0f && s_high < 1.0f) {
                float theta_high = asinf(s_high);
                float r_high = distance * tanf(2.0f * theta_high);
                r_max_sq = r_high * r_high;
            }
        }
    }

    /* Candidate pixel mask - parallelized across rows */
    int total_pixels = nx * ny;
    uint8_t *mask = (uint8_t *)calloc(total_pixels, sizeof(uint8_t));
    if (mask == NULL) {
        free(tile_mean);
        free(tile_std);
        return 0;
    }

    #pragma omp parallel for schedule(static)
    for (int y = 0; y < ny; ++y) {
        int ty = y / TILE_SIZE;
        const float *row = &data[y * nx];
        uint8_t *mask_row = &mask[y * nx];
        float ry = (y - params->beam_y) * params->pixel_size_y;
        float ry2 = ry * ry;

        for (int x = 0; x < nx; ++x) {
            float rx = (x - params->beam_x) * params->pixel_size_x;
            float r2 = rx * rx + ry2;

            if (r2 < r_min_sq || r2 > r_max_sq) {
                continue;
            }

            int tx = x / TILE_SIZE;
            int t_idx = ty * tx_count + tx;
            float v = row[x];
            float bg = tile_mean[t_idx];
            float std = tile_std[t_idx];
            float threshold = bg + params->snr_threshold * std;

            if (v > threshold && v > 0.0f) {
                mask_row[x] = 1;
            }
        }
    }

    /* Connected component analysis */
    int *queue = (int *)malloc(total_pixels * sizeof(int));
    int spot_capacity = 4096;
    int spot_count = 0;
    MxSpot *spots = (MxSpot *)malloc(spot_capacity * sizeof(MxSpot));

    if (queue == NULL || spots == NULL) {
        free(tile_mean);
        free(tile_std);
        free(mask);
        free(queue);
        free(spots);
        return 0;
    }

    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            int idx = y * nx + x;
            if (mask[idx] == 0) {
                continue;
            }

            /* Start BFS connected component */
            int head = 0;
            int tail = 0;
            queue[tail++] = idx;
            mask[idx] = 0;

            int area = 0;
            double sum_I = 0.0;
            double sum_net_I = 0.0;
            double sum_x = 0.0;
            double sum_y = 0.0;
            float max_I = -1e9f;
            float max_snr = 0.0f;

            while (head < tail) {
                int curr = queue[head++];
                int cy = curr / nx;
                int cx = curr % nx;
                float c_val = data[curr];

                int c_tile = (cy / TILE_SIZE) * tx_count + (cx / TILE_SIZE);
                float c_bg = tile_mean[c_tile];
                float c_std = tile_std[c_tile];
                float net_I = (c_val > c_bg) ? (c_val - c_bg) : 0.001f;

                area++;
                sum_I += (double)c_val;
                sum_net_I += (double)net_I;
                sum_x += (double)cx * net_I;
                sum_y += (double)cy * net_I;

                if (c_val > max_I) {
                    max_I = c_val;
                    max_snr = (c_val - c_bg) / c_std;
                }

                /* 8-connected neighbors */
                for (int dy = -1; dy <= 1; ++dy) {
                    int ny_coord = cy + dy;
                    if (ny_coord < 0 || ny_coord >= ny) continue;

                    for (int dx = -1; dx <= 1; ++dx) {
                        if (dx == 0 && dy == 0) continue;
                        int nx_coord = cx + dx;
                        if (nx_coord < 0 || nx_coord >= nx) continue;

                        int n_idx = ny_coord * nx + nx_coord;
                        if (mask[n_idx] != 0) {
                            mask[n_idx] = 0;
                            queue[tail++] = n_idx;
                        }
                    }
                }
            }

            if (area >= params->min_spot_area && area <= params->max_spot_area) {
                float cx = (sum_net_I > 0.0) ? (float)(sum_x / sum_net_I) : (float)x;
                float cy = (sum_net_I > 0.0) ? (float)(sum_y / sum_net_I) : (float)y;

                float rx = (cx - params->beam_x) * params->pixel_size_x;
                float ry = (cy - params->beam_y) * params->pixel_size_y;
                float r = sqrtf(rx * rx + ry * ry);
                float theta = 0.5f * atan2f(r, (params->distance > 0.0f ? params->distance : 100.0f));
                float sin_theta = sinf(theta);
                float d = (sin_theta > 1e-6f && params->wavelength > 0.0f)
                              ? (params->wavelength / (2.0f * sin_theta))
                              : 999.0f;

                if (params->d_min > 0.0f && d < params->d_min) {
                    continue;
                }
                if (params->d_max > 0.0f && d > params->d_max) {
                    continue;
                }

                if (spot_count >= spot_capacity) {
                    spot_capacity *= 2;
                    MxSpot *new_spots = (MxSpot *)realloc(spots, spot_capacity * sizeof(MxSpot));
                    if (new_spots == NULL) {
                        break;
                    }
                    spots = new_spots;
                }

                spots[spot_count].x = cx;
                spots[spot_count].y = cy;
                spots[spot_count].d_spacing = d;
                spots[spot_count].intensity = (float)sum_I;
                spots[spot_count].snr = max_snr;
                spot_count++;
            }
        }
    }

    /* Sort spots descending by intensity */
    if (spot_count > 1) {
        qsort(spots, spot_count, sizeof(MxSpot), compare_spots_desc);
    }

    /* Copy to out_spots if provided */
    if (out_spots != NULL && max_spots > 0) {
        int copy_count = (spot_count < max_spots) ? spot_count : max_spots;
        memcpy(out_spots, spots, copy_count * sizeof(MxSpot));
    }

    free(tile_mean);
    free(tile_std);
    free(mask);
    free(queue);
    free(spots);

    return spot_count;
}

int mxspots_score_spots(
    const MxSpot *spots,
    int spot_count,
    MxScoreResult *out_score
) {
    if (out_score == NULL) {
        return -1;
    }

    if (spots == NULL || spot_count <= 0) {
        out_score->spot_count = 0;
        out_score->avg_snr = 0.0f;
        out_score->d_min = 999.0f;
        out_score->percentage_indexed = 0.0f;
        return 0;
    }

    out_score->spot_count = spot_count;
    double sum_snr = 0.0;
    float min_d = 999.0f;

    for (int i = 0; i < spot_count; ++i) {
        sum_snr += spots[i].snr;
        if (spots[i].d_spacing > 0.0f && spots[i].d_spacing < min_d) {
            min_d = spots[i].d_spacing;
        }
    }

    out_score->avg_snr = (float)(sum_snr / spot_count);
    out_score->d_min = min_d;
    out_score->percentage_indexed = 0.0f;

    return 0;
}

int mxspots_score_frame(
    const float *data,
    int nx,
    int ny,
    const MxSpotsParams *params,
    MxScoreResult *out_score
) {
    if (data == NULL || nx <= 0 || ny <= 0 || params == NULL || out_score == NULL) {
        return -1;
    }

    int max_spots = 50000;
    MxSpot *spots = (MxSpot *)malloc(max_spots * sizeof(MxSpot));
    if (spots == NULL) {
        return -2;
    }

    int spot_count = mxspots_find_spots(data, nx, ny, params, spots, max_spots);
    int actual_count = (spot_count < max_spots) ? spot_count : max_spots;

    int ret = mxspots_score_spots(spots, actual_count, out_score);
    free(spots);
    return ret;
}

typedef struct {
    float x;
    float y;
    float z;
} Vec3;

typedef struct {
    Vec3 vec;
    float score;
    float length;
} CandidateBasis;

static int compare_candidate_desc(const void *a, const void *b) {
    const CandidateBasis *ca = (const CandidateBasis *)a;
    const CandidateBasis *cb = (const CandidateBasis *)b;
    if (cb->score > ca->score) return 1;
    if (cb->score < ca->score) return -1;
    return 0;
}

static inline float vec3_dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline Vec3 vec3_cross(Vec3 a, Vec3 b) {
    Vec3 c;
    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;
    return c;
}

static inline float vec3_norm(Vec3 a) {
    return sqrtf(vec3_dot(a, a));
}

int mxspots_index_spots(
    const MxSpot *spots,
    int spot_count,
    const MxSpotsParams *params,
    MxIndexResult *out_index
) {
    if (out_index == NULL) {
        return -1;
    }

    /* Initialize defaults */
    out_index->unit_cell[0] = 50.0f;
    out_index->unit_cell[1] = 50.0f;
    out_index->unit_cell[2] = 50.0f;
    out_index->unit_cell[3] = 90.0f;
    out_index->unit_cell[4] = 90.0f;
    out_index->unit_cell[5] = 90.0f;
    out_index->percentage_indexed = 0.0f;
    out_index->indexed_spot_count = 0;
    out_index->total_spot_count = spot_count;

    if (spots == NULL || spot_count <= 0 || params == NULL) {
        return 0;
    }

    float wavelength = (params->wavelength > 0.0f) ? params->wavelength : 1.0f;
    float distance = (params->distance > 0.0f) ? params->distance : 100.0f;
    float qx = (params->pixel_size_x > 0.0f) ? params->pixel_size_x : 0.075f;
    float qy = (params->pixel_size_y > 0.0f) ? params->pixel_size_y : 0.075f;
    float bx = params->beam_x;
    float by = params->beam_y;

    /* Compute reciprocal space vectors s_i = (s_x, s_y, s_z) for each spot */
    int n_spots = (spot_count < 1000) ? spot_count : 1000;
    Vec3 *s_vecs = (Vec3 *)malloc(n_spots * sizeof(Vec3));
    if (s_vecs == NULL) {
        return -2;
    }

    for (int i = 0; i < n_spots; ++i) {
        float px = (spots[i].x - bx) * qx;
        float py = (spots[i].y - by) * qy;
        float pz = distance;
        float R = sqrtf(px * px + py * py + pz * pz);
        s_vecs[i].x = px / (wavelength * R);
        s_vecs[i].y = py / (wavelength * R);
        s_vecs[i].z = (pz / R - 1.0f) / wavelength;
    }

    /* Generate candidate search directions */
    int max_dirs = 512;
    Vec3 *dirs = (Vec3 *)malloc(max_dirs * sizeof(Vec3));
    int n_dirs = 0;

    if (dirs == NULL) {
        free(s_vecs);
        return -2;
    }

    /* 1. Pairwise differences among top spots */
    int n_diff_spots = (n_spots < 40) ? n_spots : 40;
    for (int i = 0; i < n_diff_spots && n_dirs < max_dirs - 64; ++i) {
        for (int j = i + 1; j < n_diff_spots && n_dirs < max_dirs - 64; ++j) {
            Vec3 ds;
            ds.x = s_vecs[j].x - s_vecs[i].x;
            ds.y = s_vecs[j].y - s_vecs[i].y;
            ds.z = s_vecs[j].z - s_vecs[i].z;
            float norm = vec3_norm(ds);
            if (norm > 0.003f && norm < 0.25f) {
                dirs[n_dirs].x = ds.x / norm;
                dirs[n_dirs].y = ds.y / norm;
                dirs[n_dirs].z = ds.z / norm;
                n_dirs++;
            }
        }
    }

    /* 2. Uniform spherical spiral grid for completeness */
    int n_sphere = 64;
    for (int i = 0; i < n_sphere && n_dirs < max_dirs; ++i) {
        float z = (float)i / (float)(n_sphere - 1);
        float r = sqrtf(1.0f - z * z);
        float phi = (float)i * 2.39996322972865332f; /* golden angle */
        dirs[n_dirs].x = r * cosf(phi);
        dirs[n_dirs].y = r * sinf(phi);
        dirs[n_dirs].z = z;
        n_dirs++;
    }

    /* Evaluate 1D Fourier power spectrum along each direction in parallel */
    CandidateBasis *candidates = (CandidateBasis *)malloc(n_dirs * sizeof(CandidateBasis));
    if (candidates == NULL) {
        free(s_vecs);
        free(dirs);
        return -2;
    }

    int n_eval_spots = (n_spots < 150) ? n_spots : 150;

    #pragma omp parallel for schedule(static)
    for (int d = 0; d < n_dirs; ++d) {
        Vec3 u = dirs[d];
        float proj_local[150];
        for (int i = 0; i < n_eval_spots; ++i) {
            proj_local[i] = vec3_dot(s_vecs[i], u);
        }

        float best_power = 0.0f;
        float best_a = 50.0f;

        /* Scan trial cell period from 15.0 to 200.0 Angstroms */
        for (float a = 15.0f; a <= 200.0f; a += 0.5f) {
            double sum_cos = 0.0;
            double sum_sin = 0.0;
            double two_pi_a = 2.0 * M_PI * (double)a;

            for (int i = 0; i < n_eval_spots; ++i) {
                double phase = two_pi_a * (double)proj_local[i];
                sum_cos += cos(phase);
                sum_sin += sin(phase);
            }

            float power = (float)(sum_cos * sum_cos + sum_sin * sum_sin);
            if (power > best_power) {
                best_power = power;
                best_a = a;
            }
        }

        candidates[d].vec.x = best_a * u.x;
        candidates[d].vec.y = best_a * u.y;
        candidates[d].vec.z = best_a * u.z;
        candidates[d].score = best_power / ((float)n_eval_spots * (float)n_eval_spots);
        candidates[d].length = best_a;
    }

    free(dirs);

    /* Sort candidate basis vectors descending by Fourier power score */
    qsort(candidates, n_dirs, sizeof(CandidateBasis), compare_candidate_desc);

    /* Select 3 linearly independent direct space basis vectors */
    Vec3 a_vec = candidates[0].vec;
    Vec3 b_vec = {0.0f, 0.0f, 0.0f};
    Vec3 c_vec = {0.0f, 0.0f, 0.0f};
    int found_b = 0;
    int found_c = 0;

    float norm_a = vec3_norm(a_vec);
    if (norm_a < 1.0f) norm_a = 50.0f;

    for (int i = 1; i < n_dirs; ++i) {
        Vec3 v = candidates[i].vec;
        float norm_v = vec3_norm(v);
        if (norm_v < 1.0f) continue;

        float cos_angle = fabsf(vec3_dot(a_vec, v) / (norm_a * norm_v));
        if (cos_angle < 0.90f) { /* Angle > ~25 degrees from a_vec */
            b_vec = v;
            found_b = 1;
            break;
        }
    }

    if (!found_b) {
        b_vec.x = -a_vec.y;
        b_vec.y = a_vec.x;
        b_vec.z = 0.0f;
    }

    float norm_b = vec3_norm(b_vec);
    if (norm_b < 1.0f) norm_b = norm_a;

    Vec3 cross_ab = vec3_cross(a_vec, b_vec);
    float norm_cross = vec3_norm(cross_ab);

    for (int i = 1; i < n_dirs; ++i) {
        Vec3 v = candidates[i].vec;
        float norm_v = vec3_norm(v);
        if (norm_v < 1.0f) continue;

        if (norm_cross > 1e-4f) {
            float vol_frac = fabsf(vec3_dot(v, cross_ab)) / (norm_v * norm_cross);
            if (vol_frac > 0.20f) { /* Substantial component perpendicular to a-b plane */
                c_vec = v;
                found_c = 1;
                break;
            }
        }
    }

    if (!found_c) {
        if (norm_cross > 1e-4f) {
            float target_len = 0.5f * (norm_a + norm_b);
            c_vec.x = (cross_ab.x / norm_cross) * target_len;
            c_vec.y = (cross_ab.y / norm_cross) * target_len;
            c_vec.z = (cross_ab.z / norm_cross) * target_len;
        } else {
            c_vec.x = 0.0f;
            c_vec.y = 0.0f;
            c_vec.z = norm_a;
        }
    }

    float norm_c = vec3_norm(c_vec);
    if (norm_c < 1.0f) norm_c = norm_a;

    /* Compute unit cell lengths and angles */
    float a = norm_a;
    float b = norm_b;
    float c = norm_c;

    float cos_alpha = vec3_dot(b_vec, c_vec) / (b * c);
    float cos_beta = vec3_dot(a_vec, c_vec) / (a * c);
    float cos_gamma = vec3_dot(a_vec, b_vec) / (a * b);

    if (cos_alpha > 1.0f) cos_alpha = 1.0f;
    if (cos_alpha < -1.0f) cos_alpha = -1.0f;
    if (cos_beta > 1.0f) cos_beta = 1.0f;
    if (cos_beta < -1.0f) cos_beta = -1.0f;
    if (cos_gamma > 1.0f) cos_gamma = 1.0f;
    if (cos_gamma < -1.0f) cos_gamma = -1.0f;

    float alpha = (float)(acos(cos_alpha) * (180.0 / M_PI));
    float beta = (float)(acos(cos_beta) * (180.0 / M_PI));
    float gamma = (float)(acos(cos_gamma) * (180.0 / M_PI));

    /* Fractional indexing check across all spots */
    int indexed_count = 0;
    float tol = 0.20f;

    for (int i = 0; i < spot_count; ++i) {
        Vec3 s;
        if (i < n_spots) {
            s = s_vecs[i];
        } else {
            float px = (spots[i].x - bx) * qx;
            float py = (spots[i].y - by) * qy;
            float pz = distance;
            float R = sqrtf(px * px + py * py + pz * pz);
            s.x = px / (wavelength * R);
            s.y = py / (wavelength * R);
            s.z = (pz / R - 1.0f) / wavelength;
        }

        float h = vec3_dot(s, a_vec);
        float k = vec3_dot(s, b_vec);
        float l = vec3_dot(s, c_vec);

        float dh = fabsf(h - roundf(h));
        float dk = fabsf(k - roundf(k));
        float dl = fabsf(l - roundf(l));

        if (dh <= tol && dk <= tol && (found_c ? (dl <= tol || fabsf(l) <= 0.5f) : 1)) {
            indexed_count++;
        }
    }

    out_index->unit_cell[0] = a;
    out_index->unit_cell[1] = b;
    out_index->unit_cell[2] = c;
    out_index->unit_cell[3] = alpha;
    out_index->unit_cell[4] = beta;
    out_index->unit_cell[5] = gamma;
    out_index->indexed_spot_count = indexed_count;
    out_index->total_spot_count = spot_count;
    out_index->percentage_indexed = (spot_count > 0) ? (100.0f * (float)indexed_count / (float)spot_count) : 0.0f;

    free(s_vecs);
    free(candidates);
    return 0;
}

int mxspots_index_frame(
    const float *data,
    int nx,
    int ny,
    const MxSpotsParams *params,
    MxIndexResult *out_index
) {
    if (data == NULL || nx <= 0 || ny <= 0 || params == NULL || out_index == NULL) {
        return -1;
    }

    int max_spots = 50000;
    MxSpot *spots = (MxSpot *)malloc(max_spots * sizeof(MxSpot));
    if (spots == NULL) {
        return -2;
    }

    int spot_count = mxspots_find_spots(data, nx, ny, params, spots, max_spots);
    int actual_count = (spot_count < max_spots) ? spot_count : max_spots;

    int ret = mxspots_index_spots(spots, actual_count, params, out_index);
    free(spots);
    return ret;
}
