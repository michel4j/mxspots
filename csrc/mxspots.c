#define MXSPOTS_EXPORTS
#include "mxspots.h"

int mxspots_get_version(void) {
    return 100; /* Version 1.0.0 encoded as 100 */
}

int mxspots_ping(const MxSpotsParams *params) {
    if (params == NULL) {
        return -1;
    }
    /* Simple sanity check on parameters */
    if (params->snr_threshold <= 0.0f) {
        return -2;
    }
    return 0;
}
