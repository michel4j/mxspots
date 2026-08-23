# 0013 - Composite Quality Score in C Engine

To rapidly evaluate holistic diffraction quality for live beamline data streams, `mxspots` computes a unified 0-100 composite quality score directly within the C Spot Engine (`mxspots_score_frame`). Integrating spot count efficiency, 95th percentile resolution ($d_{95}$), Signal-to-Noise Ratio (SNR), Bragg regularity, indexing percentage, and ice contamination penalties as a native C formula minimizes Python boundary crossings and enables zero-allocation high-speed classification.
