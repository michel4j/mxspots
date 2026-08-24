#define MXSPOTS_EXPORTS
#include "mxspots.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DEFAULT_BG_HALF_WIDTH 15
#define DEFAULT_PEAK_HALF_WIDTH 3
#define DEFAULT_SPOTS_CAPACITY 50000
#define DEFAULT_CCL_CAPACITY 65536
#define MAX_RUNS_PER_ROW 4096

typedef struct {
    int parent;
    int area;
    double sum_I;
    double sum_net_I;
    double sum_x;
    double sum_y;
    float max_I;
    float max_snr;
} CclComponent;

typedef struct {
    int x_start;
    int x_end;
    int label;
} Run;

struct MxSpotsContext {
    int max_nx;
    int max_ny;
    double *sat_sum;
    double *sat_sum_sq;
    int *sat_count;
    int8_t *mask;
    MxSpot *spots_buf;
    int spots_capacity;
    CclComponent *ccl_nodes;
    int ccl_capacity;
};

typedef struct {
    double *sum;
    double *sum_sq;
    int *count;
    int nx;
    int ny;
    int stride;
} IntegralImage;

static inline void sat_query_box(
    const IntegralImage *sat,
    int x0, int y0, int x1, int y1,
    double *out_sum, double *out_sum_sq, int *out_count
) {
    int stride = sat->stride;
    int i_br = y1 * stride + x1;
    int i_tr = y0 * stride + x1;
    int i_bl = y1 * stride + x0;
    int i_tl = y0 * stride + x0;

    *out_sum = sat->sum[i_br] - sat->sum[i_tr] - sat->sum[i_bl] + sat->sum[i_tl];
    *out_sum_sq = sat->sum_sq[i_br] - sat->sum_sq[i_tr] - sat->sum_sq[i_bl] + sat->sum_sq[i_tl];
    *out_count = sat->count[i_br] - sat->count[i_tr] - sat->count[i_bl] + sat->count[i_tl];
}

static inline void sat_query_annulus(
    const IntegralImage *sat,
    int x, int y,
    int r_out, int r_in,
    float *out_mean, float *out_std
) {
    int x0_out = (x - r_out > 0) ? (x - r_out) : 0;
    int y0_out = (y - r_out > 0) ? (y - r_out) : 0;
    int x1_out = (x + r_out + 1 < sat->nx) ? (x + r_out + 1) : sat->nx;
    int y1_out = (y + r_out + 1 < sat->ny) ? (y + r_out + 1) : sat->ny;

    int x0_in = (x - r_in > 0) ? (x - r_in) : 0;
    int y0_in = (y - r_in > 0) ? (y - r_in) : 0;
    int x1_in = (x + r_in + 1 < sat->nx) ? (x + r_in + 1) : sat->nx;
    int y1_in = (y + r_in + 1 < sat->ny) ? (y + r_in + 1) : sat->ny;

    double s_out, sq_out;
    int n_out;
    sat_query_box(sat, x0_out, y0_out, x1_out, y1_out, &s_out, &sq_out, &n_out);

    double s_in, sq_in;
    int n_in;
    sat_query_box(sat, x0_in, y0_in, x1_in, y1_in, &s_in, &sq_in, &n_in);

    double s_ann = s_out - s_in;
    double sq_ann = sq_out - sq_in;
    int n_ann = n_out - n_in;

    if (n_ann > 5) {
        double mean = s_ann / n_ann;
        double var = (sq_ann / n_ann) - (mean * mean);
        double std = (var > 0.0) ? sqrt(var) : 1.0;
        *out_mean = (float)mean;
        *out_std = (std < 1.0) ? 1.0f : (float)std;
    } else if (n_out > 0) {
        double mean = s_out / n_out;
        double var = (sq_out / n_out) - (mean * mean);
        double std = (var > 0.0) ? sqrt(var) : 1.0;
        *out_mean = (float)mean;
        *out_std = (std < 1.0) ? 1.0f : (float)std;
    } else {
        *out_mean = 0.0f;
        *out_std = 1.0f;
    }
}

static inline int ccl_find_root(CclComponent *nodes, int label) {
    int root = label;
    while (nodes[root].parent != root) {
        root = nodes[root].parent;
    }
    int curr = label;
    while (curr != root) {
        int next = nodes[curr].parent;
        nodes[curr].parent = root;
        curr = next;
    }
    return root;
}

static inline void ccl_union(CclComponent *nodes, int root1, int root2) {
    if (root1 == root2) return;
    nodes[root2].parent = root1;
    nodes[root1].area += nodes[root2].area;
    nodes[root1].sum_I += nodes[root2].sum_I;
    nodes[root1].sum_net_I += nodes[root2].sum_net_I;
    nodes[root1].sum_x += nodes[root2].sum_x;
    nodes[root1].sum_y += nodes[root2].sum_y;
    if (nodes[root2].max_I > nodes[root1].max_I) {
        nodes[root1].max_I = nodes[root2].max_I;
        nodes[root1].max_snr = nodes[root2].max_snr;
    }
}

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

static int compare_float_asc(const void *a, const void *b) {
    float fa = *(const float *)a;
    float fb = *(const float *)b;
    if (fa < fb) return -1;
    if (fa > fb) return 1;
    return 0;
}

typedef struct {
    float x;
    float y;
    float z;
} Vec3;

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

typedef struct {
    Vec3 vec;
    int count;
} DiffCluster;

static int compare_diff_cluster_desc(const void *a, const void *b) {
    const DiffCluster *ca = (const DiffCluster *)a;
    const DiffCluster *cb = (const DiffCluster *)b;
    if (cb->count > ca->count) return 1;
    if (cb->count < ca->count) return -1;
    return 0;
}

static inline int uf_find(int *parent, int i) {
    int root = i;
    while (parent[root] != root) root = parent[root];
    int curr = i;
    while (curr != root) {
        int next = parent[curr];
        parent[curr] = root;
        curr = next;
    }
    return root;
}

static inline void uf_union(int *parent, int *size, int i, int j) {
    int root1 = uf_find(parent, i);
    int root2 = uf_find(parent, j);
    if (root1 != root2) {
        if (size[root1] < size[root2]) {
            parent[root1] = root2;
            size[root2] += size[root1];
        } else {
            parent[root2] = root1;
            size[root1] += size[root2];
        }
    }
}

int mxspots_analyze_regularity(
    const MxSpot *spots,
    int spot_count,
    const MxSpotsParams *params,
    float *out_bragg_percent,
    int *out_bragg_spots,
    float *out_avg_intensity,
    int *out_num_lattices,
    float *out_d_min
) {
    if (out_bragg_percent != NULL) *out_bragg_percent = 0.0f;
    if (out_bragg_spots != NULL) *out_bragg_spots = 0;
    if (out_avg_intensity != NULL) *out_avg_intensity = 0.0f;
    if (out_num_lattices != NULL) *out_num_lattices = 0;
    if (out_d_min != NULL) *out_d_min = 999.0f;

    if (spots == NULL || spot_count < 5 || params == NULL) {
        return 0;
    }

    float wavelength = (params->wavelength > 0.0f) ? params->wavelength : 1.0f;
    float distance = (params->distance > 0.0f) ? params->distance : 200.0f;
    float qx = (params->pixel_size_x > 0.0f) ? params->pixel_size_x : 0.075f;
    float qy = (params->pixel_size_y > 0.0f) ? params->pixel_size_y : 0.075f;
    float bx = params->beam_x;
    float by = params->beam_y;

    int n_spots = (spot_count < 2000) ? spot_count : 2000;
    Vec3 *s_vecs = (Vec3 *)malloc(n_spots * sizeof(Vec3));
    if (s_vecs == NULL) {
        return -1;
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

    /* 1. Extract pairwise difference vectors from top spots */
    int n_diff_spots = (n_spots < 150) ? n_spots : 150;
    int max_diffs = (n_diff_spots * (n_diff_spots - 1)) / 2;
    Vec3 *diffs = (Vec3 *)malloc(max_diffs * sizeof(Vec3));
    if (diffs == NULL) {
        free(s_vecs);
        return -1;
    }

    int num_diffs = 0;
    for (int i = 0; i < n_diff_spots; ++i) {
        for (int j = i + 1; j < n_diff_spots; ++j) {
            Vec3 ds;
            ds.x = s_vecs[j].x - s_vecs[i].x;
            ds.y = s_vecs[j].y - s_vecs[i].y;
            ds.z = s_vecs[j].z - s_vecs[i].z;
            float norm = vec3_norm(ds);
            /* Corresponds to real space repeats between ~10 A and ~250 A */
            if (norm >= 0.004f && norm <= 0.085f) {
                if (ds.x < 0.0f || (ds.x == 0.0f && ds.y < 0.0f) || (ds.x == 0.0f && ds.y == 0.0f && ds.z < 0.0f)) {
                    ds.x = -ds.x;
                    ds.y = -ds.y;
                    ds.z = -ds.z;
                }
                diffs[num_diffs++] = ds;
            }
        }
    }

    if (num_diffs < 3) {
        free(s_vecs);
        free(diffs);
        return 0;
    }

    /* 2. Count recurrence / density for difference vectors */
    int *recurrence = (int *)calloc(num_diffs, sizeof(int));
    if (recurrence == NULL) {
        free(s_vecs);
        free(diffs);
        return -1;
    }

    float tol_diff = 0.0045f;
    for (int i = 0; i < num_diffs; ++i) {
        for (int j = 0; j < num_diffs; ++j) {
            float dx = diffs[i].x - diffs[j].x;
            float dy = diffs[i].y - diffs[j].y;
            float dz = diffs[i].z - diffs[j].z;
            if (sqrtf(dx*dx + dy*dy + dz*dz) <= tol_diff) {
                recurrence[i]++;
            }
        }
    }

    /* 3. Extract dominant recurring difference vectors */
    DiffCluster *clusters = (DiffCluster *)malloc(num_diffs * sizeof(DiffCluster));
    if (clusters == NULL) {
        free(s_vecs);
        free(diffs);
        free(recurrence);
        return -1;
    }

    for (int i = 0; i < num_diffs; ++i) {
        clusters[i].vec = diffs[i];
        clusters[i].count = recurrence[i];
    }
    free(diffs);
    free(recurrence);

    qsort(clusters, num_diffs, sizeof(DiffCluster), compare_diff_cluster_desc);

    Vec3 test_vectors[128];
    int num_test_vectors = 0;
    int max_test_vectors = 64;
    int min_recurrence = (n_spots < 150) ? 2 : 3;

    /* Collect top distinct recurring difference vectors */
    for (int i = 0; i < num_diffs && num_test_vectors < max_test_vectors; ++i) {
        if (clusters[i].count < min_recurrence) break;

        Vec3 cand = clusters[i].vec;
        float cand_norm = vec3_norm(cand);
        if (cand_norm < 1e-4f) continue;

        int duplicate = 0;
        for (int b = 0; b < num_test_vectors; ++b) {
            float dx = cand.x - test_vectors[b].x;
            float dy = cand.y - test_vectors[b].y;
            float dz = cand.z - test_vectors[b].z;
            if (sqrtf(dx*dx + dy*dy + dz*dz) < 0.0035f) {
                duplicate = 1;
                break;
            }
        }

        if (!duplicate) {
            test_vectors[num_test_vectors++] = cand;
        }
    }

    free(clusters);

    if (num_test_vectors == 0) {
        free(s_vecs);
        return 0;
    }

    /* 4. Build spot recurrence adjacency graph using Union-Find */
    int *uf_parent = (int *)malloc(n_spots * sizeof(int));
    int *uf_size = (int *)malloc(n_spots * sizeof(int));
    if (uf_parent == NULL || uf_size == NULL) {
        free(s_vecs);
        free(uf_parent);
        free(uf_size);
        return -1;
    }

    for (int i = 0; i < n_spots; ++i) {
        uf_parent[i] = i;
        uf_size[i] = 1;
    }

    float tol_edge = 0.0050f;
    for (int i = 0; i < n_spots; ++i) {
        for (int j = i + 1; j < n_spots; ++j) {
            Vec3 ds;
            ds.x = s_vecs[j].x - s_vecs[i].x;
            ds.y = s_vecs[j].y - s_vecs[i].y;
            ds.z = s_vecs[j].z - s_vecs[i].z;

            int is_edge = 0;
            for (int b = 0; b < num_test_vectors; ++b) {
                Vec3 u = test_vectors[b];
                /* Check +u and -u */
                float ex1 = ds.x - u.x;
                float ey1 = ds.y - u.y;
                float ez1 = ds.z - u.z;
                if (sqrtf(ex1*ex1 + ey1*ey1 + ez1*ez1) <= tol_edge) {
                    is_edge = 1;
                    break;
                }

                float ex2 = ds.x + u.x;
                float ey2 = ds.y + u.y;
                float ez2 = ds.z + u.z;
                if (sqrtf(ex2*ex2 + ey2*ey2 + ez2*ez2) <= tol_edge) {
                    is_edge = 1;
                    break;
                }
            }

            if (is_edge) {
                uf_union(uf_parent, uf_size, i, j);
            }
        }
    }

    /* 5. Connected component extraction and lattice quantification */
    int bragg_spots = 0;
    int num_lattices = 0;
    int min_comp_size = 5;
    double sum_bragg_intensity = 0.0;

    float *bragg_d = (float *)malloc(n_spots * sizeof(float));
    int bragg_d_count = 0;

    for (int i = 0; i < n_spots; ++i) {
        if (uf_parent[i] == i) {
            if (uf_size[i] >= min_comp_size) {
                num_lattices++;
            }
        }
        int root = uf_find(uf_parent, i);
        if (uf_size[root] >= min_comp_size) {
            bragg_spots++;
            sum_bragg_intensity += (double)spots[i].intensity;
            if (spots[i].d_spacing > 0.0f && spots[i].d_spacing < 900.0f && bragg_d != NULL) {
                bragg_d[bragg_d_count++] = spots[i].d_spacing;
            }
        }
    }

    float bragg_percent = (spot_count > 0) ? (100.0f * (float)bragg_spots / (float)spot_count) : 0.0f;
    if (bragg_percent > 100.0f) bragg_percent = 100.0f;

    float avg_intensity = (bragg_spots > 0) ? (float)(sum_bragg_intensity / (double)bragg_spots) : 0.0f;

    /* Compute 98th percentile resolution limit exclusively from confirmed Bragg spots */
    float bragg_d_min = 999.0f;
    if (bragg_d_count > 0 && bragg_d != NULL) {
        qsort(bragg_d, bragg_d_count, sizeof(float), compare_float_asc);
        int k = (int)(0.02f * (float)(bragg_d_count - 1));
        if (k < 0) k = 0;
        if (k >= bragg_d_count) k = bragg_d_count - 1;
        bragg_d_min = bragg_d[k];
    }
    if (bragg_d != NULL) {
        free(bragg_d);
    }

    if (out_bragg_percent != NULL) *out_bragg_percent = bragg_percent;
    if (out_bragg_spots != NULL) *out_bragg_spots = bragg_spots;
    if (out_avg_intensity != NULL) *out_avg_intensity = avg_intensity;
    if (out_num_lattices != NULL) *out_num_lattices = num_lattices;
    if (out_d_min != NULL) *out_d_min = bragg_d_min;

    free(s_vecs);
    free(uf_parent);
    free(uf_size);
    return 0;
}

MxSpotsContext *mxspots_create_context(int max_nx, int max_ny) {
    if (max_nx <= 0 || max_ny <= 0) {
        return NULL;
    }

    MxSpotsContext *ctx = (MxSpotsContext *)calloc(1, sizeof(MxSpotsContext));
    if (ctx == NULL) {
        return NULL;
    }

    ctx->max_nx = max_nx;
    ctx->max_ny = max_ny;
    ctx->spots_capacity = DEFAULT_SPOTS_CAPACITY;
    ctx->ccl_capacity = DEFAULT_CCL_CAPACITY;

    int stride = max_nx + 1;
    size_t sat_entries = (size_t)(max_ny + 1) * stride;
    size_t total_pixels = (size_t)max_nx * max_ny;

    ctx->sat_sum = (double *)calloc(sat_entries, sizeof(double));
    ctx->sat_sum_sq = (double *)calloc(sat_entries, sizeof(double));
    ctx->sat_count = (int *)calloc(sat_entries, sizeof(int));
    ctx->mask = (int8_t *)calloc(total_pixels, sizeof(int8_t));
    ctx->spots_buf = (MxSpot *)malloc(ctx->spots_capacity * sizeof(MxSpot));
    ctx->ccl_nodes = (CclComponent *)malloc(ctx->ccl_capacity * sizeof(CclComponent));

    if (ctx->sat_sum == NULL || ctx->sat_sum_sq == NULL || ctx->sat_count == NULL ||
        ctx->mask == NULL || ctx->spots_buf == NULL || ctx->ccl_nodes == NULL) {
        mxspots_free_context(ctx);
        return NULL;
    }

    return ctx;
}

void mxspots_free_context(MxSpotsContext *ctx) {
    if (ctx == NULL) {
        return;
    }
    free(ctx->sat_sum);
    free(ctx->sat_sum_sq);
    free(ctx->sat_count);
    free(ctx->mask);
    free(ctx->spots_buf);
    free(ctx->ccl_nodes);
    free(ctx);
}

int mxspots_find_spots_ctx(
    MxSpotsContext *ctx,
    const float *data,
    int nx,
    int ny,
    const MxSpotsParams *params,
    MxSpot *out_spots,
    int max_spots
) {
    if (ctx == NULL || data == NULL || nx <= 0 || ny <= 0 || params == NULL) {
        return 0;
    }

    if (nx > ctx->max_nx || ny > ctx->max_ny) {
        return 0;
    }

    int stride = nx + 1;
    IntegralImage sat;
    sat.sum = ctx->sat_sum;
    sat.sum_sq = ctx->sat_sum_sq;
    sat.count = ctx->sat_count;
    sat.nx = nx;
    sat.ny = ny;
    sat.stride = stride;

    /* Pass 1: Compute cumulative row integrals in parallel */
    #pragma omp parallel for schedule(static)
    for (int y = 0; y < ny; ++y) {
        double r_sum = 0.0;
        double r_sum_sq = 0.0;
        int r_count = 0;
        int row_offset = (y + 1) * stride;
        const float *row = &data[y * nx];

        for (int x = 0; x < nx; ++x) {
            float v = row[x];
            if (v >= 0.0f && !isnan(v) && !isinf(v)) {
                r_sum += (double)v;
                r_sum_sq += (double)v * v;
                r_count += 1;
            }
            ctx->sat_sum[row_offset + (x + 1)] = r_sum;
            ctx->sat_sum_sq[row_offset + (x + 1)] = r_sum_sq;
            ctx->sat_count[row_offset + (x + 1)] = r_count;
        }
    }

    /* Pass 2: Compute cumulative column integrals in parallel */
    #pragma omp parallel for schedule(static)
    for (int x = 1; x <= nx; ++x) {
        double c_sum = 0.0;
        double c_sum_sq = 0.0;
        int c_count = 0;
        for (int y = 1; y <= ny; ++y) {
            int idx = y * stride + x;
            c_sum += ctx->sat_sum[idx];
            c_sum_sq += ctx->sat_sum_sq[idx];
            c_count += ctx->sat_count[idx];
            ctx->sat_sum[idx] = c_sum;
            ctx->sat_sum_sq[idx] = c_sum_sq;
            ctx->sat_count[idx] = c_count;
        }
    }

    /* Safe geometry parameter fallback defaults */
    float distance = (params->distance > 0.0f) ? params->distance : 200.0f;
    float wavelength = (params->wavelength > 0.0f) ? params->wavelength : 1.0f;
    float pixel_size_x = (params->pixel_size_x > 0.0f) ? params->pixel_size_x : 0.075f;
    float pixel_size_y = (params->pixel_size_y > 0.0f) ? params->pixel_size_y : 0.075f;
    float beam_x = (params->beam_x != 0.0f) ? params->beam_x : 0.0f;
    float beam_y = (params->beam_y != 0.0f) ? params->beam_y : 0.0f;

    /* Compute resolution radial limits */
    float r_min_sq = 0.0f;
    float r_max_sq = 1e30f;

    if (wavelength > 0.0f && distance > 0.0f) {
        if (params->d_max > 0.0f) {
            float s_low = wavelength / (2.0f * params->d_max);\
            if (s_low > 0.0f && s_low < 1.0f) {
                float theta_low = asinf(s_low);\
                float r_low = distance * tanf(2.0f * theta_low);\
                r_min_sq = r_low * r_low;\
            }
        }
        if (params->d_min > 0.0f) {
            float s_high = wavelength / (2.0f * params->d_min);\
            if (s_high > 0.0f && s_high < 1.0f) {
                float theta_high = asinf(s_high);\
                float r_high = distance * tanf(2.0f * theta_high);\
                r_max_sq = r_high * r_high;\
            }
        }
    }

    /* Clear and populate candidate pixel mask */
    int total_pixels = nx * ny;
    int8_t *mask = ctx->mask;
    memset(mask, 0, total_pixels * sizeof(int8_t));

    #pragma omp parallel for schedule(static)
    for (int y = 0; y < ny; ++y) {
        const float *row = &data[y * nx];
        int8_t *mask_row = &mask[y * nx];
        float ry = (y - beam_y) * pixel_size_y;
        float ry2 = ry * ry;

        for (int x = 0; x < nx; ++x) {
            float rx = (x - beam_x) * pixel_size_x;
            float r2 = rx * rx + ry2;

            if (r2 < r_min_sq || r2 > r_max_sq) {
                continue;
            }

            /* Check masked ice rings */
            int is_masked = 0;
            for (int k = 0; k < params->num_masked_rings; ++k) {
                if (r2 >= params->masked_rings_r2[k][0] && r2 <= params->masked_rings_r2[k][1]) {
                    is_masked = 1;
                    break;
                }
            }
            if (is_masked) {
                continue;
            }

            float v = row[x];
            if (v <= 0.0f) {
                continue;
            }

            float bg, std;
            sat_query_annulus(&sat, x, y, DEFAULT_BG_HALF_WIDTH, DEFAULT_PEAK_HALF_WIDTH, &bg, &std);

            float threshold = bg + params->snr_threshold * std;
            if (v > threshold) {
                mask_row[x] = 1; /* Candidate strong pixel */
            }
        }
    }

    /* Streaming Two-Pass Run-Length Connected Component Labeling */
    CclComponent *nodes = ctx->ccl_nodes;
    int num_nodes = 0;
    int max_nodes = ctx->ccl_capacity;

    Run prev_runs[MAX_RUNS_PER_ROW];
    Run curr_runs[MAX_RUNS_PER_ROW];
    int num_prev_runs = 0;

    for (int y = 0; y < ny; ++y) {
        const float *row = &data[y * nx];
        const int8_t *mask_row = &mask[y * nx];
        int num_curr_runs = 0;
        int x = 0;

        while (x < nx) {
            if (mask_row[x] == 1) {
                int x_start = x;
                while (x < nx && mask_row[x] == 1) {
                    x++;
                }
                int x_end = x - 1;

                if (num_curr_runs < MAX_RUNS_PER_ROW && num_nodes < max_nodes) {
                    int run_area = 0;
                    double run_sum_I = 0.0;
                    double run_sum_net_I = 0.0;
                    double run_sum_x = 0.0;
                    double run_sum_y = 0.0;
                    float run_max_I = -1e9f;
                    float run_max_snr = 0.0f;

                    for (int rx = x_start; rx <= x_end; ++rx) {
                        float v = row[rx];
                        float bg, std;
                        sat_query_annulus(&sat, rx, y, DEFAULT_BG_HALF_WIDTH, DEFAULT_PEAK_HALF_WIDTH, &bg, &std);
                        float net_I = (v > bg) ? (v - bg) : 0.001f;
                        float snr = (v - bg) / std;

                        run_area++;
                        run_sum_I += (double)v;
                        run_sum_net_I += (double)net_I;
                        run_sum_x += (double)rx * net_I;
                        run_sum_y += (double)y * net_I;

                        if (v > run_max_I) {
                            run_max_I = v;
                            run_max_snr = snr;
                        }
                    }

                    int matched_root = -1;
                    for (int p = 0; p < num_prev_runs; ++p) {
                        /* 8-connectivity with previous row: x_start <= prev.x_end + 1 && x_end >= prev.x_start - 1 */
                        if (x_start <= prev_runs[p].x_end + 1 && x_end >= prev_runs[p].x_start - 1) {
                            int prev_root = ccl_find_root(nodes, prev_runs[p].label);
                            if (matched_root == -1) {
                                matched_root = prev_root;
                            } else if (matched_root != prev_root) {
                                ccl_union(nodes, matched_root, prev_root);
                            }
                        }
                    }

                    int label;
                    if (matched_root == -1) {
                        label = num_nodes++;
                        nodes[label].parent = label;
                        nodes[label].area = run_area;
                        nodes[label].sum_I = run_sum_I;
                        nodes[label].sum_net_I = run_sum_net_I;
                        nodes[label].sum_x = run_sum_x;
                        nodes[label].sum_y = run_sum_y;
                        nodes[label].max_I = run_max_I;
                        nodes[label].max_snr = run_max_snr;
                    } else {
                        label = matched_root;
                        nodes[label].area += run_area;
                        nodes[label].sum_I += run_sum_I;
                        nodes[label].sum_net_I += run_sum_net_I;
                        nodes[label].sum_x += run_sum_x;
                        nodes[label].sum_y += run_sum_y;
                        if (run_max_I > nodes[label].max_I) {
                            nodes[label].max_I = run_max_I;
                            nodes[label].max_snr = run_max_snr;
                        }
                    }

                    curr_runs[num_curr_runs].x_start = x_start;
                    curr_runs[num_curr_runs].x_end = x_end;
                    curr_runs[num_curr_runs].label = label;
                    num_curr_runs++;
                }
            } else {
                x++;
            }
        }

        if (num_curr_runs > 0) {
            memcpy(prev_runs, curr_runs, num_curr_runs * sizeof(Run));
            num_prev_runs = num_curr_runs;
        } else {
            num_prev_runs = 0;
        }
    }

    /* Second Pass: Filter and extract valid root components */
    int spot_count = 0;
    MxSpot *spots = ctx->spots_buf;
    int spot_capacity = ctx->spots_capacity;

    for (int i = 0; i < num_nodes; ++i) {
        if (nodes[i].parent == i) { /* Root node */
            int area = nodes[i].area;
            if (area >= params->min_spot_area && area <= params->max_spot_area) {
                double net_I = nodes[i].sum_net_I;
                float cx = (net_I > 0.0) ? (float)(nodes[i].sum_x / net_I) : 0.0f;
                float cy = (net_I > 0.0) ? (float)(nodes[i].sum_y / net_I) : 0.0f;

                float rx = (cx - beam_x) * pixel_size_x;
                float ry = (cy - beam_y) * pixel_size_y;
                float r2 = rx * rx + ry * ry;

                int spot_masked = 0;
                for (int k = 0; k < params->num_masked_rings; ++k) {
                    if (r2 >= params->masked_rings_r2[k][0] && r2 <= params->masked_rings_r2[k][1]) {
                        spot_masked = 1;
                        break;
                    }
                }
                if (spot_masked) {
                    continue;
                }

                float r = sqrtf(r2);
                float theta = 0.5f * atan2f(r, distance);
                float sin_theta = sinf(theta);
                float d = (sin_theta > 1e-6f && wavelength > 0.0f)
                              ? (wavelength / (2.0f * sin_theta))
                              : 999.0f;

                if (params->d_min > 0.0f && d < params->d_min) {
                    continue;
                }
                if (params->d_max > 0.0f && d > params->d_max) {
                    continue;
                }

                if (spot_count < spot_capacity) {
                    spots[spot_count].x = cx;
                    spots[spot_count].y = cy;
                    spots[spot_count].d_spacing = d;
                    spots[spot_count].intensity = (float)nodes[i].sum_I;
                    spots[spot_count].snr = nodes[i].max_snr;
                    spot_count++;
                }
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

    return spot_count;
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

    MxSpotsContext *ctx = mxspots_create_context(nx, ny);
    if (ctx == NULL) {
        return 0;
    }

    int spot_count = mxspots_find_spots_ctx(ctx, data, nx, ny, params, out_spots, max_spots);
    mxspots_free_context(ctx);
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
        out_score->bragg_spots = 0;
        out_score->bragg_percent = 0.0f;
        out_score->avg_intensity = 0.0f;
        out_score->avg_snr = 0.0f;
        out_score->d_min = 999.0f;
        out_score->ice_score = 0.0f;
        out_score->num_ice_rings = 0;
        out_score->num_lattices = 0;
        out_score->score = 0.0f;
        return 0;
    }

    out_score->spot_count = spot_count;
    double sum_snr = 0.0;
    for (int i = 0; i < spot_count; ++i) {
        sum_snr += spots[i].snr;
    }
    out_score->avg_snr = (float)(sum_snr / spot_count);

    int n_bragg = out_score->bragg_spots;
    if (n_bragg <= 0 && out_score->bragg_percent > 0.0f) {
        n_bragg = (int)roundf((out_score->bragg_percent / 100.0f) * (float)spot_count);
    }

    /* Strict gate: if no Bragg spots, resolution is unresolved (999.0) and score is 0.0 */
    if (n_bragg <= 0) {
        out_score->d_min = 999.0f;
        out_score->score = 0.0f;
        return 0;
    }

    /* If d_min was not already provided from regularity analysis (e.g. out_score->d_min <= 0.0f or >= 900.0f), compute from spots */
    if (out_score->d_min <= 0.0f || out_score->d_min >= 900.0f) {
        float *d_spacings = (float *)malloc(spot_count * sizeof(float));
        int valid_d_count = 0;
        for (int i = 0; i < spot_count; ++i) {
            if (spots[i].d_spacing > 0.0f && spots[i].d_spacing < 900.0f) {
                if (d_spacings != NULL) {
                    d_spacings[valid_d_count++] = spots[i].d_spacing;
                }
            }
        }
        if (valid_d_count > 0 && d_spacings != NULL) {
            qsort(d_spacings, valid_d_count, sizeof(float), compare_float_asc);
            int k = (int)(0.02f * (float)(valid_d_count - 1));
            if (k < 0) k = 0;
            if (k >= valid_d_count) k = valid_d_count - 1;
            out_score->d_min = d_spacings[k];
        } else {
            out_score->d_min = 999.0f;
        }
        if (d_spacings != NULL) {
            free(d_spacings);
        }
    }

    float d_98 = out_score->d_min;

    float p_bragg = (out_score->bragg_percent > 0.0f) ? out_score->bragg_percent : (100.0f * (float)n_bragg / (float)spot_count);
    if (p_bragg > 100.0f) p_bragg = 100.0f;
    if (p_bragg < 0.0f) p_bragg = 0.0f;

    float snr = (out_score->avg_snr > 0.0f) ? out_score->avg_snr : 0.0f;

    float s_res = 0.0f;
    if (d_98 < 4.0f) {
        s_res = (4.0f - d_98) / (4.0f - 1.2f);
        if (s_res > 1.0f) s_res = 1.0f;
        if (s_res < 0.0f) s_res = 0.0f;
    }

    float p_ice = 0.8f * (float)out_score->num_ice_rings;
    if (out_score->ice_score > 2.0f) {
        p_ice += 0.4f * (out_score->ice_score - 2.0f);
    }

    float z = -5.5f + 0.85f * logf(1.0f + (float)n_bragg)
                    + 1.20f * (p_bragg / 100.0f)
                    + 0.50f * logf(1.0f + snr)
                    + 0.60f * s_res
                    - p_ice;

    float score = 100.0f / (1.0f + expf(-z));
    if (score > 100.0f) score = 100.0f;
    if (score < 0.0f) score = 0.0f;

    out_score->score = score;
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

    /* 1. Ice ring detection */
    MxIceResult ice_res;
    ice_res.num_rings = 0;
    ice_res.ice_score = 0.0f;
    mxspots_detect_ice(data, nx, ny, params, &ice_res);
    out_score->ice_score = ice_res.ice_score;
    out_score->num_ice_rings = ice_res.num_rings;

    /* 2. Spot finding */
    int max_spots = 50000;
    MxSpot *spots = (MxSpot *)malloc(max_spots * sizeof(MxSpot));
    if (spots == NULL) {
        return -2;
    }

    int spot_count = mxspots_find_spots(data, nx, ny, params, spots, max_spots);
    int actual_count = (spot_count < max_spots) ? spot_count : max_spots;

    /* 3. Regularity analysis */
    out_score->bragg_spots = 0;
    out_score->bragg_percent = 0.0f;
    out_score->avg_intensity = 0.0f;
    out_score->num_lattices = 0;
    out_score->d_min = 999.0f;
    if (actual_count >= 5) {
        mxspots_analyze_regularity(
            spots,
            actual_count,
            params,
            &out_score->bragg_percent,
            &out_score->bragg_spots,
            &out_score->avg_intensity,
            &out_score->num_lattices,
            &out_score->d_min
        );
    }

    int ret = mxspots_score_spots(spots, actual_count, out_score);
    free(spots);
    return ret;
}

/* Canonical Ice Ih d-spacings in Angstroms */
static const float CANONICAL_ICE_D_SPACINGS[] = {
    3.897f, /* (100) */
    3.669f, /* (002) */
    3.441f, /* (101) */
    2.249f, /* (102) */
    2.071f, /* (110) */
    1.918f  /* (103) */
};

int mxspots_detect_ice(
    const float *data,
    int nx,
    int ny,
    const MxSpotsParams *params,
    MxIceResult *out_result
) {
    if (data == NULL || nx <= 0 || ny <= 0 || params == NULL || out_result == NULL) {
        return -1;
    }

    out_result->num_rings = 0;
    out_result->ice_score = 0.0f;

    float distance = (params->distance > 0.0f) ? params->distance : 200.0f;
    float wavelength = (params->wavelength > 0.0f) ? params->wavelength : 1.0f;
    float pixel_size_x = (params->pixel_size_x > 0.0f) ? params->pixel_size_x : 0.075f;
    float pixel_size_y = (params->pixel_size_y > 0.0f) ? params->pixel_size_y : 0.075f;
    float beam_x = (params->beam_x != 0.0f) ? params->beam_x : (nx / 2.0f);
    float beam_y = (params->beam_y != 0.0f) ? params->beam_y : (ny / 2.0f);
    float sensitivity = (params->snr_threshold > 0.0f) ? params->snr_threshold : 3.0f;

    float q_avg = 0.5f * (pixel_size_x + pixel_size_y);
    if (q_avg <= 0.0f) q_avg = 0.075f;

    /* Compute maximum pixel radius */
    float dx0 = beam_x * pixel_size_x;
    float dx1 = ((float)nx - 1.0f - beam_x) * pixel_size_x;
    float dy0 = beam_y * pixel_size_y;
    float dy1 = ((float)ny - 1.0f - beam_y) * pixel_size_y;
    float max_dx = (dx0 > dx1) ? dx0 : dx1;
    float max_dy = (dy0 > dy1) ? dy0 : dy1;
    float max_r_mm = sqrtf(max_dx * max_dx + max_dy * max_dy);
    float max_r_px = max_r_mm / q_avg;
    int num_bins = (int)ceilf(max_r_px) + 2;
    if (num_bins < 100) num_bins = 100;

    double *radial_sum = (double *)calloc(num_bins, sizeof(double));
    int *radial_count = (int *)calloc(num_bins, sizeof(int));
    if (radial_sum == NULL || radial_count == NULL) {
        free(radial_sum);
        free(radial_count);
        return -2;
    }

    #pragma omp parallel
    {
        double *t_sum = (double *)calloc(num_bins, sizeof(double));
        int *t_cnt = (int *)calloc(num_bins, sizeof(int));

        if (t_sum != NULL && t_cnt != NULL) {
            #pragma omp for schedule(static)
            for (int y = 0; y < ny; ++y) {
                const float *row = &data[y * nx];
                float ry = (y - beam_y) * pixel_size_y;
                float ry2 = ry * ry;

                for (int x = 0; x < nx; ++x) {
                    float v = row[x];
                    if (v <= 0.0f || isnan(v) || isinf(v)) continue;

                    float rx = (x - beam_x) * pixel_size_x;
                    float r_mm = sqrtf(rx * rx + ry2);
                    int bin = (int)roundf(r_mm / q_avg);
                    if (bin >= 0 && bin < num_bins) {
                        t_sum[bin] += (double)v;
                        t_cnt[bin] += 1;
                    }
                }
            }

            #pragma omp critical
            {
                for (int b = 0; b < num_bins; ++b) {
                    radial_sum[b] += t_sum[b];
                    radial_count[b] += t_cnt[b];
                }
            }
        }

        free(t_sum);
        free(t_cnt);
    }

    float *radial_mean = (float *)malloc(num_bins * sizeof(float));
    if (radial_mean == NULL) {
        free(radial_sum);
        free(radial_count);
        return -2;
    }

    for (int b = 0; b < num_bins; ++b) {
        radial_mean[b] = (radial_count[b] > 5) ? (float)(radial_sum[b] / radial_count[b]) : 0.0f;
    }

    int num_canonical = sizeof(CANONICAL_ICE_D_SPACINGS) / sizeof(CANONICAL_ICE_D_SPACINGS[0]);
    float max_overall_score = 0.0f;

    for (int i = 0; i < num_canonical && out_result->num_rings < MXSPOTS_MAX_MASKED_RINGS; ++i) {
        float d_ice = CANONICAL_ICE_D_SPACINGS[i];
        float sin_theta = wavelength / (2.0f * d_ice);
        if (sin_theta <= 0.0f || sin_theta >= 1.0f) continue;

        float theta = asinf(sin_theta);
        float r_mm = distance * tanf(2.0f * theta);
        float r_px = r_mm / q_avg;
        int b_center = (int)roundf(r_px);

        if (b_center < 10 || b_center >= num_bins - 10) continue;

        /* Calculate local background in annulus around b_center */
        double bg_sum = 0.0;
        double bg_sum_sq = 0.0;
        int bg_cnt = 0;

        for (int b = b_center - 15; b <= b_center + 15; ++b) {
            if (b < 0 || b >= num_bins) continue;
            if (b >= b_center - 4 && b <= b_center + 4) continue;
            if (radial_count[b] > 5) {
                float val = radial_mean[b];
                bg_sum += val;
                bg_sum_sq += val * val;
                bg_cnt++;
            }
        }

        if (bg_cnt < 5) continue;

        double bg_mean = bg_sum / bg_cnt;
        double var = (bg_sum_sq / bg_cnt) - (bg_mean * bg_mean);
        double bg_std = (var > 0.0) ? sqrt(var) : 1.0;
        if (bg_std < 0.1) bg_std = 1.0;

        /* Find peak in range [b_center - 3, b_center + 3] */
        float peak_val = 0.0f;
        for (int b = b_center - 3; b <= b_center + 3; ++b) {
            if (b >= 0 && b < num_bins) {
                if (radial_mean[b] > peak_val) {
                    peak_val = radial_mean[b];
                }
            }
        }

        float snr = (float)((peak_val - bg_mean) / bg_std);
        if (snr > max_overall_score) {
            max_overall_score = snr;
        }

        if (snr >= sensitivity) {
            float delta_r_px = 5.0f;
            float r_low_mm = (r_px - delta_r_px > 1.0f) ? ((r_px - delta_r_px) * q_avg) : (1.0f * q_avg);
            float r_high_mm = (r_px + delta_r_px) * q_avg;

            float theta_low = 0.5f * atan2f(r_low_mm, distance);
            float sin_t_low = sinf(theta_low);
            float d_max_ring = (sin_t_low > 1e-6f) ? (wavelength / (2.0f * sin_t_low)) : 30.0f;

            float theta_high = 0.5f * atan2f(r_high_mm, distance);
            float sin_t_high = sinf(theta_high);
            float d_min_ring = (sin_t_high > 1e-6f) ? (wavelength / (2.0f * sin_t_high)) : 1.0f;

            int ring_idx = out_result->num_rings;
            out_result->rings[ring_idx].d_spacing = d_ice;
            out_result->rings[ring_idx].d_min = d_min_ring;
            out_result->rings[ring_idx].d_max = d_max_ring;
            out_result->rings[ring_idx].score = snr;
            out_result->num_rings++;
        }
    }

    out_result->ice_score = max_overall_score;

    free(radial_sum);
    free(radial_count);
    free(radial_mean);

    return 0;
}
