#define MXSPOTS_EXPORTS
#include "mxspots.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define TILE_SIZE 32

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

    /* Pass 1: Compute tile background mean and std */
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

    /* Candidate pixel mask */
    int total_pixels = nx * ny;
    uint8_t *mask = (uint8_t *)calloc(total_pixels, sizeof(uint8_t));
    if (mask == NULL) {
        free(tile_mean);
        free(tile_std);
        return 0;
    }

    for (int y = 0; y < ny; ++y) {
        int ty = y / TILE_SIZE;
        const float *row = &data[y * nx];
        uint8_t *mask_row = &mask[y * nx];

        for (int x = 0; x < nx; ++x) {
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
