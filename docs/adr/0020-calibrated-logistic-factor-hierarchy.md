# 0020 - Calibrated Logistic Factor Hierarchy and Quality Scoring

## Context
ADR 0016 introduced a Hybrid Gated-Logistic composite scoring function with a bounded post-logistic ice penalty. However, in the initial calibration ($w_N = 0.85$, $w_S = 0.50$, $w_P = 1.20$, $w_{\text{res}} = 0.60$), Bragg spot count ($N_B$) accounted for ~56% of total logit variation while Signal-to-Noise Ratio (SNR) accounted for only ~21%. This caused noisy frames with large numbers of weak spots to receive high scores, while strongly diffracting frames with fewer reflections were under-scored. Furthermore, the Bragg percentage coefficient ($w_P = 1.20$) nominally exceeded the spot count and SNR weights.

A principled recalibration establishes a clear three-tier factor hierarchy:
1. **Tier 1 (Major Factors, ~82% dynamic share)**: Bragg spot count ($N_B$) and Average SNR ($\text{SNR}$) co-equally drive the score.
2. **Tier 2 (Secondary Factor, ~12% dynamic share)**: Bragg percentage ($P_B$) modulates the score based on crystalline diffraction purity.
3. **Tier 3 (Tertiary Factor, ~6% dynamic share)**: 98th percentile resolution limit ($d_{98}$) provides refinement for high-resolution diffraction.

## Decision
Adopt a recalibrated logit formulation in the C Spot Engine:

1. **Hard Gate**: If confirmed Bragg spot count $N_B = 0$, $\text{Score} = 0.0$.
2. **Base Logistic Score**: If $N_B > 0$,
   $$\text{Score}_{\text{raw}} = \frac{100}{1 + \exp(-z)}$$
   $$\text{Score} = \text{clamp}\left(\text{Score}_{\text{raw}} - P_{\text{ice}}, 0.0, 100.0\right)$$
3. **Calibrated Logit Formula**:
   $$z = w_0 + w_N \ln(1 + N_B) + w_S \ln(1 + \text{SNR}) + w_P (P_B / 100) + w_{\text{res}} s_{\text{res}}$$
   where calibrated coefficients are:
   - $w_0 = -6.50$ (calibrated baseline intercept)
   - $w_N = 0.75$ (Bragg spots count weight, ~44.2% dynamic share)
   - $w_S = 1.20$ (Signal-to-Noise Ratio weight, ~38.2% dynamic share)
   - $w_P = 1.40$ (Bragg percentage weight, ~11.9% dynamic share)
   - $w_{\text{res}} = 0.50$ (Resolution limit weight, ~5.6% dynamic share)
   - $s_{\text{res}} = \text{clamp}\left(\frac{4.0 - d_{98}}{4.0 - 1.2}, 0.0, 1.0\right)$
4. **Post-Logistic Ice Penalty ($P_{\text{ice}} \in [0.0, 10.0]$)**:
   $$P_{\text{ice}} = \text{clamp}\left(2.0 \times N_{\text{ice}} + 1.0 \times \max(0.0, I_{\text{ice}} - 2.0), 0.0, 10.0\right)$$

## Consequences
- Noisy diffraction frames with high spot counts but low SNR (e.g. SNR $\approx 3.5$) are properly suppressed in quality score.
- High-intensity, well-diffracting frames with fewer spots are awarded higher scores due to strong SNR weighting.
- Full backwards compatibility is preserved for downstream Python APIs and C ABI.
