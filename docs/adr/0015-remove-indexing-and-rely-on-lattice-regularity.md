# 0015 - Remove Indexing and Rely on Lattice Regularity for Bragg Distinction

## Context
Reciprocal space FFT-based indexing was originally included to compute indexed spot percentages and unit cell parameters. However, full indexing adds algorithmic complexity and is sensitive to incomplete single-frame coverage. Lattice regularity analysis via difference-vector recurrence graphs and connected components already provides robust, fast, and rotation-invariant identification of periodic Bragg reflections without assuming or determining unit cell geometry.

## Decision
Completely remove reciprocal lattice indexing (`indexer.py`, `mxspots_index_spots`, `mxspots_index_frame`, `MxIndexResult`, `IndexResult`, and `mxspots.index` CLI). Rely solely on difference-vector recurrence graph connected components to distinguish Bragg reflections from non-Bragg spots. In `ScoreResult` and `MxScoreResult`, replace indexing and old regularity fields with `bragg_spots`, `bragg_percent`, and `avg_intensity` (average intensity of Bragg spots).
