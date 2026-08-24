# 0016 - Hybrid Gated-Logistic Quality Scoring with Resolution Limit

## Context
The previous linear scoring formula suffered from hard boundary clipping and could assign non-zero scores to noisy frames with no periodic crystalline diffraction. A smooth logistic formulation with a strict gate on Bragg reflection presence ensures that only genuine crystalline diffraction patterns receive positive quality scores.

## Decision
Implement a Hybrid Gated-Logistic composite scoring function in the C Spot Engine:
1. **Hard Gate**: If Bragg spot count $N_B = 0$, $\text{Score} = 0.0$.
2. **Logistic Curve**: If $N_B > 0$, $\text{Score} = \frac{100}{1 + \exp(-z)}$.
3. **Logit Formula**:
   $$z = w_0 + w_N \ln(1 + N_B) + w_P (P_B / 100) + w_S \ln(1 + \text{SNR}) + w_{\text{res}} s_{\text{res}} - P_{\text{ice}}$$
   where calibrated weights are:
   - $w_0 = -5.5$ (base intercept)
   - $w_N = 0.85$ (Bragg spots count weight)
   - $w_P = 1.20$ (Bragg percentage weight)
   - $w_S = 0.50$ (Signal-to-noise ratio weight)
   - $w_{\text{res}} = 0.60$ (Resolution limit weight)
   - $s_{\text{res}} = \text{clamp}\left(\frac{4.0 - d_{98}}{4.0 - 1.2}, 0.0, 1.0\right)$
   - $P_{\text{ice}} = 0.8 \times N_{\text{ice}} + 0.4 \times \max(0.0, I_{\text{ice}} - 2.0)$
