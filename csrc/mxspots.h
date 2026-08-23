#ifndef MXSPOTS_H
#define MXSPOTS_H

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32) || defined(__CYGWIN__)
  #ifdef MXSPOTS_EXPORTS
    #define MXSPOTS_API __declspec(dllexport)
  #else
    #define MXSPOTS_API __declspec(dllimport)
  #endif
#else
  #if __GNUC__ >= 4
    #define MXSPOTS_API __attribute__ ((visibility ("default")))
  #else
    #define MXSPOTS_API
  #endif
#endif

#include <stdint.h>

#define MXSPOTS_MAX_MASKED_RINGS 16

/**
 * Parameters controlling spot detection, resolution range, and detector geometry.
 */
typedef struct {
    float snr_threshold;       /* Signal-to-noise ratio threshold */
    int min_spot_area;         /* Minimum connected pixels */
    int max_spot_area;         /* Maximum connected pixels */
    float beam_x;              /* Beam center X in pixels */
    float beam_y;              /* Beam center Y in pixels */
    float pixel_size_x;        /* Detector pixel size X in mm */
    float pixel_size_y;        /* Detector pixel size Y in mm */
    float distance;            /* Sample-to-detector distance in mm */
    float wavelength;          /* Incident X-ray wavelength in Angstroms */
    float d_min;               /* High-resolution limit in Angstroms (0 = unbounded) */
    float d_max;               /* Low-resolution limit in Angstroms */
    int num_masked_rings;      /* Number of active masked radial/resolution rings */
    float masked_rings_r2[MXSPOTS_MAX_MASKED_RINGS][2]; /* [min_r2, max_r2] in mm^2 */
} MxSpotsParams;

/**
 * Detected spot representation.
 */
typedef struct {
    float x;                   /* Centroid X in pixel coordinates */
    float y;                   /* Centroid Y in pixel coordinates */
    float d_spacing;           /* Bragg resolution d-spacing in Angstroms */
    float intensity;           /* Integrated spot intensity */
    float snr;                 /* Peak SNR */
} MxSpot;

/**
 * Diffraction frame quality score metrics.
 */
typedef struct {
    int spot_count;            /* Number of detected spots */
    float avg_snr;             /* Average signal-to-noise ratio */
    float d_min;               /* Maximum observable resolution in Angstroms */
    float percentage_indexed;  /* Indexable fraction percentage (0.0 if not indexed) */
} MxScoreResult;

/**
 * Reciprocal lattice indexing result.
 */
typedef struct {
    float unit_cell[6];        /* a, b, c (Angstroms), alpha, beta, gamma (degrees) */
    float percentage_indexed;  /* Fraction of spots indexed to candidate lattice (0-100%) */
    int indexed_spot_count;    /* Number of spots conforming to lattice */
    int total_spot_count;      /* Total spots evaluated */
} MxIndexResult;

/**
 * Reusable execution scratch context for zero-allocation batch spot finding.
 */
typedef struct MxSpotsContext MxSpotsContext;

/**
 * Returns the library version code.
 */
MXSPOTS_API int mxspots_get_version(void);

/**
 * Validates parameter configuration.
 */
MXSPOTS_API int mxspots_ping(const MxSpotsParams *params);

/**
 * Allocates a reusable context with scratch buffers for frames up to max_nx x max_ny.
 */
MXSPOTS_API MxSpotsContext *mxspots_create_context(int max_nx, int max_ny);

/**
 * Frees a previously created context and its internal scratch buffers.
 */
MXSPOTS_API void mxspots_free_context(MxSpotsContext *ctx);

/**
 * Core spot finding routine using a pre-allocated scratch context (zero allocations).
 */
MXSPOTS_API int mxspots_find_spots_ctx(
    MxSpotsContext *ctx,
    const float *data,
    int nx,
    int ny,
    const MxSpotsParams *params,
    MxSpot *out_spots,
    int max_spots
);

/**
 * Core spot finding routine. Detects diffraction spots in a 2D float image frame.
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
 * Computes diffraction frame quality score from a pre-calculated spot list.
 */
MXSPOTS_API int mxspots_score_spots(
    const MxSpot *spots,
    int spot_count,
    MxScoreResult *out_score
);

/**
 * End-to-end frame scoring routine: finds spots and computes quality metrics.
 */
MXSPOTS_API int mxspots_score_frame(
    const float *data,
    int nx,
    int ny,
    const MxSpotsParams *params,
    MxScoreResult *out_score
);

/**
 * Indexes reciprocal lattice directly from a detected spot list.
 */
MXSPOTS_API int mxspots_index_spots(
    const MxSpot *spots,
    int spot_count,
    const MxSpotsParams *params,
    MxIndexResult *out_index
);

/**
 * Indexes reciprocal lattice directly from a raw 2D diffraction image frame.
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
