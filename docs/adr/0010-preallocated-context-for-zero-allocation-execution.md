# 0010 - Pre-allocated Context for Zero-Allocation Frame Ingestion

To support high-throughput beamline data streams without per-frame dynamic memory allocation overhead, `mxspots` provides `MxSpotsContext` / `SpotFinderContext`. Reusable Summed-Area Tables, candidate pixel masks, and component buffers are allocated once and reused across continuous frame acquisitions.
