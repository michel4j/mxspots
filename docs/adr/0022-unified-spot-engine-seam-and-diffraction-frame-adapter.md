# 0022 - Unified Spot Engine Seam and DiffractionFrame Adapter

## Context
Prior to this architectural deepening, several core pipelines exhibited leaky seams and shallow interfaces:
1. **Spot Scoring Coordination**: Python `score_spots` was forced to coordinate two separate C entry points (`mxspots_analyze_regularity` and `mxspots_score_spots`), manually marshalling seven ctypes pointer arguments across the foreign function interface.
2. **Coupled Frame Ingestion**: Frame loading, format detection (mxio, YAML, synthetic, duck typing), and detector geometry merging were implemented in a loose procedural helper (`extract_frame_and_params`) in `spotfinder.py`. `scorer.py` depended directly on `spotfinder.py` for this helper and for ice detection, coupling two independent workflows.
3. **Split Ice Ring Masking**: Ice detection and post-mask validation loops were fragmented across Python and C. `score_data` ran Python ice detection while `mxspots_score_frame` repeated ice detection in C. `findspots_data` ran a secondary Python spot filtering loop over all detected spots.

## Decision

1. **Deepened Spot Engine Scoring Seam**:
   - `mxspots_score_spots` in C now accepts `(const MxSpot *spots, int spot_count, const MxSpotsParams *params, MxScoreResult *out_score)`.
   - When `params` is provided and `out_score->bragg_spots <= 0`, `mxspots_score_spots` internally executes reciprocal lattice regularity analysis and computes composite quality scoring behind a single interface.
   - Preserves backward compatibility: if `params == NULL` or Bragg metrics are already populated, existing metrics are retained.

2. **DiffractionFrame Adapter**:
   - Introduced `DiffractionFrame` adapter in `models.py` with `DiffractionFrame.from_source(image_source, params=None)`.
   - Encapsulates 2D NumPy array validation, format-specific ingestion (`mxio.DataSet`, `mxio.ImageFrame`, `SyntheticFrame`, YAML specs), detector geometry resolution, and frame indexing.
   - `findspots`, `score`, and `detect_ice_rings` consume `DiffractionFrame`, decoupling `scorer.py` from `spotfinder.py` internals.
   - `extract_frame_and_params` is retained as a thin compatibility wrapper.

3. **Unified C Spot Engine Ice Ring Masking**:
   - `MxScoreResult` embeds `MxIceRing rings[MXSPOTS_MAX_MASKED_RINGS]`. `mxspots_score_frame` returns complete ice ring metadata in a single pass without redundant Python detection.
   - `mxspots_find_spots_ctx` and `mxspots_find_spots` accept `MxIceResult *out_ice` and an `auto_ice_mask` parameter flag.
   - 1D radial ice profiling, empirical annulus expansion, radial interval consolidation, CCL masking, and post-mask spot exclusion reside entirely within the C Spot Engine. Python-level spot iteration loops are eliminated.

## Consequences
- **Locality**: Format ingestion is isolated in `DiffractionFrame`; regularity analysis and ice masking policies are concentrated in the Spot Engine.
- **Leverage**: Callers invoke a single entry point for spot finding and scoring without intermediate pointer marshalling or multi-stage parameter slicing.
- **Testability**: Tests exercise the full scoring and ice-exclusion pipeline through deep, consolidated interfaces.
