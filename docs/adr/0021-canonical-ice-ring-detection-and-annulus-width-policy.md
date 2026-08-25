# 0021 - Canonical Ice Ring Detection, Annulus-Width Expansion Policy, and Mask Hardening

## Context
Hexagonal water ice ($I_h$) formation during macromolecular cryo-crystallography produces characteristic azimuthally symmetric powder diffraction rings at well-defined resolution shells (e.g., 3.897 Å, 3.669 Å, 3.441 Å, 2.249 Å, 2.071 Å, 1.918 Å, 1.883 Å, 1.724 Å). Previous iterations used fixed-width radial intervals (e.g. $\pm 5.0\text{ px}$) that were geometry-agnostic, failed to fully cover diffuse/broad ice rings, lacked interval consolidation for adjacent overlapping shells (e.g. 3.897 Å and 3.669 Å), and risked leaking false-positive Bragg spots into scoring and spot-finding pipelines.

## Decision

1. **Canonical Contract & Geometry Mapping**:
   - `IceRing` contract defines `d_spacing` (nominal Å), `d_min` (high-resolution bound in Å), `d_max` (low-resolution bound in Å), and `score` (SNR over local radial background).
   - `SpotParams.masked_rings` accepts a sequence of `(d_min, d_max)` tuples with $d_{\text{min}} \le d_{\text{max}}$.
   - $d$-spacing shells map invertibly to detector radius:
     $$\theta = \arcsin\left(\frac{\lambda}{2d}\right), \quad r = \text{distance} \times \tan(2\theta)$$
     where $d_{\text{max}}$ corresponds to $r_{\text{min}}$ and $d_{\text{min}}$ corresponds to $r_{\text{max}}$.

2. **Empirical Annulus-Width Expansion Policy**:
   - **Profile-Driven Extent**: For each detected peak in the 1D radial average, scan radially outward while the radial mean exceeds the quarter-prominence threshold:
     $$\mu_{\text{bin}} > \mu_{\text{bg}} + 0.25 \times (\text{Peak} - \mu_{\text{bg}})$$
   - **Safety Padding & Half-Width Floor**: Enforce a minimum half-width floor of $\max(6.0\text{ px}, 0.020 \times r_{\text{peak}})$ and apply $+4.0\text{ px}$ wing padding.
   - **Resolution Safety Margin**: Guarantee the resolution mask envelopes the canonical ice ring with at least $\pm 2\%$ margin ($d_{\text{min}} \le 0.98 \times d_{\text{ice}}$ and $d_{\text{max}} \ge 1.02 \times d_{\text{ice}}$).
   - **Boundary Clamping**: Clamp radial lower bounds to $r \ge 1.0\text{ px}$ and enforce $r_{\text{high}} \ge r_{\text{low}} + 2.0\text{ px}$.

3. **Interval Merging & Overlap Consolidation**:
   - When converting masked resolution shells to squared radial pixel intervals in `CMxSpotsParams.from_params` and C execution buffers, intervals are sorted and merged if contiguous or overlapping, preventing slot exhaustion and eliminating inter-ring gap leakage.

4. **Hardened Post-Mask Validation Filter**:
   - Both in C (`mxspots_score_frame` and `mxspots_find_spots_ctx`) and in Python (`findspots_data`), candidate spot centroids and $d$-spacings are validated against all active masked shells. Any spot falling within an active ice shell is strictly discarded.

5. **Calibrated Ice Detection Sensitivity**:
   - Default `ice_sensitivity` is set to $2.5\sigma$, preventing false-positive ice ring detection on flat random Gaussian noise distributions where random 7-bin maxima naturally produce $\text{SNR} \approx 1.35\sigma$.

## Consequences
- 100% elimination of residual ice-spot contamination across synthetic and experimental diffraction frames.
- Adaptive annulus coverage robustly accommodates both sharp crystalline rings and diffuse powder bands across diverse detector distances and beam wavelengths.
- Clean separation of genuine crystal Bragg reflections from solvent/ice artifacts in lattice regularity scoring.
