# 0016 - Hybrid Gated-Logistic Quality Scoring with Resolution Limit

> _Note: The logit coefficients in this ADR were recalibrated in [ADR 0020](0020-calibrated-logistic-factor-hierarchy.md) to establish Bragg spot count and SNR as the primary co-equal drivers of composite quality scoring._

## Context
The previous linear scoring formula suffered from hard boundary clipping and could assign non-zero scores to noisy frames with no periodic crystalline diffraction. A smooth logistic formulation with a strict gate on Bragg reflection presence ensures that only genuine crystalline diffraction patterns receive positive quality scores. Furthermore, placing ice penalties directly inside the exponential logit $z$ caused disproportionately punitive drops (>30-40 points) on typical frames.

## Decision
Implement a Hybrid Gated-Logistic composite scoring function with a post-logistic bounded ice penalty in the C Spot Engine:
1. **Hard Gate**: If Bragg spot count $N_B = 0$, $\text{Score} = 0.0$.
2. **Base Logistic Score**: If $N_B > 0$,
   $$\text{Score}_{\text{raw}} = \frac{100}{1 + \exp(-z)}$$
   $$\text{Score} = \text{clamp}\left(\text{Score}_{\text{raw}} - P_{\text{ice}}, 0.0, 100.0\right)$$
3. **Logit Formula**:
   $$z = w_0 + w_N \ln(1 + N_B) + w_P (P_B / 100) + w_S \ln(1 + \text{SNR}) + w_{\text{res}} s_{\text{res}}$$
   where calibrated weights are:
   - $w_0 = -5.5$ (base intercept)
   - $w_N = 0.85$ (Bragg spots count weight)
   - $w_P = 1.20$ (Bragg percentage weight)
   - $w_S = 0.50$ (Signal-to-noise ratio weight)
   - $w_{\text{res}} = 0.60$ (Resolution limit weight)
   - $s_{\text{res}} = \text{clamp}\left(\frac{4.0 - d_{98}}{4.0 - 1.2}, 0.0, 1.0\right)$
4. **Post-Logistic Ice Penalty ($P_{\text{ice}} \in [0.0, 10.0]$)**:
   $$P_{\text{ice}} = \text{clamp}\left(2.0 \times N_{\text{ice}} + 1.0 \times \max(0.0, I_{\text{ice}} - 2.0),\, 0.0,\, 10.0\right)$$
   This ensures that the presence of maximum ice contamination penalizes the frame quality score by at most 10 points (10% of maximum score).
