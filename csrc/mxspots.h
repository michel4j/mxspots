#ifndef MXSPOTS_H
#define MXSPOTS_H

#include <stdint.h>
#include <stddef.h>

#ifdef _WIN32
    #ifdef MXSPOTS_EXPORTS
        #define MXSPOTS_API __declspec(dllexport)
    #else
        #define MXSPOTS_API __declspec(dllimport)
    #endif
#else
    #define MXSPOTS_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float snr_threshold;       /* Signal-to-noise ratio threshold */
    int min_spot_area;         /* Minimum connected pixels for a spot */
    int max_spot_area;         /* Maximum connected pixels for a spot */
    float beam_x;              /* Beam center X (pixels) */
    float beam_y;              /* Beam center Y (pixels) */
    float pixel_size_x;        /* Pixel size X (mm) */
    float pixel_size_y;        /* Pixel size Y (mm) */
    float distance;            /* Detector distance (mm) */
    float wavelength;          /* X-ray wavelength (Angstroms) */
    float d_min;               /* High-resolution cutoff (Angstroms, 0 for unbounded) */
    float d_max;               /* Low-resolution cutoff (Angstroms, e.g. 30.0) */
} MxSpotsParams;

typedef struct {
    float x;
    float y;
    float d_spacing;
    float intensity;
    float snr;
} MxSpot;

typedef struct {
    int spot_count;
    float avg_snr;
    float d_min;
    float percentage_indexed;
} MxScoreResult;

typedef struct {
    float unit_cell[6];        /* a, b, c (Angstroms), alpha, beta, gamma (degrees) */
    float percentage_indexed;  /* 0.0 - 100.0 % */
    int indexed_spot_count;    /* Number of spots matching lattice */
    int total_spot_count;      /* Total number of spots evaluated */
} MxIndexResult;

/**
 * Returns the version number of libmxspots.
 */
MXSPOTS_API int mxspots_get_version(void);

/**
 * Simple ping function validating parameter struct transfer.
 * Returns 0 on success, non-zero if params pointer is NULL.
 */
MXSPOTS_API int mxspots_ping(const MxSpotsParams *params);

/**
 * Detect spots in a 2D float32 diffraction image.
 *
 * @param data       Pointer to row-major 2D float32 image array (size nx * ny).
 * @param nx         Image width in pixels.
 * @param ny         Image height in pixels.
 * @param params     Spot finding parameters and experimental geometry.
 * @param out_spots  Buffer to store detected spots (can be NULL to only count).
 * @param max_spots  Maximum capacity of out_spots.
 * @return           Total number of detected spots.
 */
MXSPOTS_API int mxspots_find_spots(
    const float *data,
    int nx,
    int ny,
    const MxSpotsParams *params,
    MxSpot *out_spots,
    int max_spots
);

/**
 * Compute quality score metrics from a detected spot list.
 *
 * @param spots      Array of detected spots.
 * @param spot_count Number of spots in array.
 * @param out_score  Pointer to MxScoreResult struct to receive quality metrics.
 * @return           0 on success, non-zero on error.
 */
MXSPOTS_API int mxspots_score_spots(
    const MxSpot *spots,
    int spot_count,
    MxScoreResult *out_score
);

/**
 * Compute quality score metrics directly from a 2D float32 diffraction image.
 *
 * @param data       Pointer to row-major 2D float32 image array (size nx * ny).
 * @param nx         Image width in pixels.
 * @param ny         Image height in pixels.
 * @param params     Spot finding parameters and experimental geometry.
 * @param out_score  Pointer to MxScoreResult struct to receive quality metrics.
 * @return           0 on success, non-zero on error.
 */
MXSPOTS_API int mxspots_score_frame(
    const float *data,
    int nx,
    int ny,
    const MxSpotsParams *params,
    MxScoreResult *out_score
);

/**
 * Index a list of detected spots using 1D/3D reciprocal lattice FFT search.
 *
 * @param spots      Array of detected spots.
 * @param spot_count Number of spots in array.
 * @param params     Experimental geometry parameters.
 * @param out_index  Pointer to MxIndexResult struct to receive unit cell and percentage indexed.
 * @return           0 on success, non-zero on error.
 */
MXSPOTS_API int mxspots_index_spots(
    const MxSpot *spots,
    int spot_count,
    const MxSpotsParams *params,
    MxIndexResult *out_index
);

/**
 * Find spots and index lattice directly from a 2D float32 diffraction image.
 *
 * @param data       Pointer to row-major 2D float32 image array (size nx * ny).
 * @param nx         Image width in pixels.
 * @param ny         Image height in pixels.
 * @param params     Spot finding parameters and experimental geometry.
 * @param out_index  Pointer to MxIndexResult struct to receive unit cell and percentage indexed.
 * @return           0 on success, non-zero on error.
 */
MXSPOTS_API int mxspots_index_frame(
    const float *data,
    int nx,
    int ny,
    const MxSpotsParams *params,
    MxIndexResult *out_index
);

#ifdef __cplusplus
}
#endif

#endif /* MXSPOTS_H */
