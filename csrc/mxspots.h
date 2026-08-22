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

/**
 * Returns the version number of libmxspots.
 */
MXSPOTS_API int mxspots_get_version(void);

/**
 * Simple ping function validating parameter struct transfer.
 * Returns 0 on success, non-zero if params pointer is NULL.
 */
MXSPOTS_API int mxspots_ping(const MxSpotsParams *params);

#ifdef __cplusplus
}
#endif

#endif /* MXSPOTS_H */
