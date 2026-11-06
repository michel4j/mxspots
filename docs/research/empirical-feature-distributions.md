# Empirical Feature Distributions Across Experimental Datasets

**Ticket**: michel4j/mxspots#89  
**Parent**: [Wayfinder Map: Diffraction Frame Quality Score Dynamic Range Discrimination (michel4j/mxspots#88)](https://github.com/michel4j/mxspots/issues/88)  
**Branch**: `research/empirical-feature-distributions`  
**Date**: November 2026  

---

## 1. Executive Summary

This study establishes the empirical distributions of spot finding, regularity metrics, resolution limits, ice contamination, and composite score logit factors across real-world macromolecular crystallography (MX) detector frames. Analysis covers 181 experimental frames across three distinct crystal systems:
- **`lyso1`** (Hen egg-white lysozyme, Pilatus 6M detector, 123 frames across 3 raster sweeps)
- **`insul-2`** (Bovine insulin, high-resolution cryo sweep, 8 frames with varying ice contamination)
- **`thaumatin`** (Thaumatococcus daniellii thaumatin, 50 frames sampled evenly across 500 frames)

### Key Empirical Findings

1. **Diffracting Frame Feature Medians vs 95th Percentiles**:
   - Confirmed Bragg spots ($N_B$): Median = 150.50 spots; 95th percentile = 627.25 spots (max = 815.00 spots). For the `lyso1` sweep alone, median $N_B = 318.00$ and 95th percentile = 692.25 spots.
   - Bragg-isolated SNR ($\text{SNR}$): Median = 64.52; 95th percentile = 768.08 (max = 2,266.78). In `lyso1`, SNR regularly exceeds 200–800 on core crystal passes.
   - Bragg purity ($P_B$): Median = 88.64%; 95th percentile = 98.15% (mean = 88.05%).
   - Resolution limit ($d_{98}$): Median = 1.23 Å; 95th percentile = 2.10 Å (range: 1.09 Å to 2.41 Å).

2. **Severe Score Saturation on Diffracting Frames**:
   - The current logistic scoring function ($z = -6.50 + 0.75 \ln(1+N_B) + 1.20 \ln(1+\text{SNR}) + 1.40 (P_B/100) + 0.50 s_{\text{res}}$) produces logit index values $z \in [3.56, 7.50]$ for $>75\%$ of diffracting frames (median $z = 3.89$ aggregate; median $z = 6.54$ in `lyso1`).
   - Because $\sigma(z) = \frac{100}{1 + e^{-z}}$ asymptotes toward 100 for $z > 3.0$, raw quality scores collapse into an extremely narrow band of **97.00 to 99.94** across almost all diffracting crystal frames, entirely eliminating discrimination between mediocre and stellar diffraction.

3. **Diagnosis of Blank Frame $\text{SNR} \ge 66,000$ Artifact**:
   - **Root Cause**: An unmasked persistent hot pixel cluster at $(y=1135, x=2219\text{--}2220)$ on the detector exhibits raw intensities of $\sim 73,800$ and $\sim 133,300$. On blank frames, the surrounding annular background has zero intensity and zero variance, causing local background standard deviation to default to $\text{std} = 1.0\text{f}$. This produces a raw spot SNR of $(133,327 - 0) / 1.0 = 133,327$.
   - **Engine Behavior**: On diffracting frames, lattice graph recurrence filtering ($N \ge 5$) correctly discards the hot pixel as non-crystalline. On blank frames ($N_B = 0$), `csrc/mxspots.c` lines 1028–1034 fall back to calculating the average SNR across all candidate spots ($133,327 / 2 \approx 66,667$).
   - While `csrc/mxspots.c` line 1040 gates the final quality score to `0.0f`, `out_score->avg_snr` leaks the $66,000+$ value into `ScoreResult`, contaminating metadata and threatening un-gated models.

---

## 2. Experimental Datasets Overview

| Dataset | Crystal / Target | Detector Format | Total Frames Analyzed | Diffracting Frames ($N_B > 0$) | Blank Frames ($N_B = 0$) | Key Characteristics |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`lyso1`** | Hen Egg-White Lysozyme | $2463 \times 2527$ CBF | 123 | 48 (39.0%) | 75 (61.0%) | Mesh loop raster grid scans; contains crystal passes and blank air/loop passes; persistent hot pixel at $(1135, 2220)$. |
| **`insul-2`** | Bovine Insulin | $2463 \times 2527$ CBF | 8 | 8 (100.0%) | 0 (0.0%) | Ultra-high resolution ($d_{98} \approx 1.11$ Å); prominent ice ring contamination on 3/8 frames ($I_{\text{ice}}$ up to 52.2). |
| **`thaumatin`** | Thaumatin | $2463 \times 2527$ CBF | 50 (sampled / 500) | 50 (100.0%) | 0 (0.0%) | Highly uniform, well-centered rotation sweep; $d_{98} \approx 1.21$ Å; minimal ice ($I_{\text{ice}} \approx 0.16$). |
| **All Combined** | Multi-system | CBF | **181** | **106** | **75** | Comprehensive benchmarking baseline for scoring model calibration. |

---

## 3. Empirical Feature Distributions

Distributions are reported with: Minimum, 25th Percentile ($P_{25}$), Median ($P_{50}$), 75th Percentile ($P_{75}$), 95th Percentile ($P_{95}$), Maximum, and Mean $\pm$ Standard Deviation.

### 3.1 `lyso1` Specific Sweep: `lyso_1-0616-094123_*.cbf`

This 49-frame sweep represents a complete raster grid scan across a crystal mounted on a mesh loop. Frames 0001–0007, 0014–0015, 0021–0022, 0028–0030, 0035–0037, and 0040–0049 capture blank air/loop, while intervening intervals capture crystal passes.

#### Diffracting Frames ($N_B > 0$, $N = 22$)
| Metric | Min | 25th | Median | 75th | 95th | Max | Mean ± Std |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Total Spots ($N_{\text{spots}}$) | 70.00 | 162.75 | 327.00 | 608.75 | 706.05 | 827.00 | 368.50 ± 234.60 |
| Bragg Spots ($N_B$) | 66.00 | 154.75 | 318.00 | 597.50 | 692.25 | 815.00 | 357.59 ± 233.46 |
| Average SNR ($\text{SNR}$) | 20.15 | 60.22 | 141.27 | 252.52 | 417.00 | 806.68 | 192.36 ± 181.47 |
| Bragg Purity ($P_B$, %) | 86.24 | 93.94 | 96.49 | 98.10 | 98.55 | 98.55 | 95.59 ± 3.01 |
| Resolution Limit ($d_{98}$, Å) | 1.91 | 1.92 | 1.98 | 2.08 | 2.25 | 2.41 | 2.03 ± 0.13 |
| Ice Profile Peak ($I_{\text{ice}}$) | 0.02 | 0.03 | 0.04 | 0.18 | 0.44 | 0.53 | 0.12 ± 0.14 |
| Detected Ice Rings ($N_{\text{ice}}$) | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 ± 0.00 |
| Ice Penalty ($P_{\text{ice}}$) | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 ± 0.00 |
| Total Quality Score | 87.52 | 97.87 | 99.76 | 99.87 | 99.89 | 99.92 | 98.06 ± 3.21 |

#### Blank Frames ($N_B = 0$, $N = 27$)
| Metric | Min | 25th | Median | 75th | 95th | Max | Mean ± Std |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Total Spots ($N_{\text{spots}}$) | 2.00 | 2.00 | 2.00 | 2.00 | 2.70 | 3.00 | 2.07 ± 0.26 |
| Bragg Spots ($N_B$) | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 ± 0.00 |
| Average SNR ($\text{SNR}$) | 44,506.38 | 66,684.84 | 66,780.39 | 66,829.51 | 66,854.14 | 66,858.69 | 65,124.76 ± 5821.93 |
| Bragg Purity ($P_B$, %) | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 ± 0.00 |
| Resolution Limit ($d_{98}$, Å) | 999.00 | 999.00 | 999.00 | 999.00 | 999.00 | 999.00 | 999.00 ± 0.00 |
| Ice Profile Peak ($I_{\text{ice}}$) | 0.01 | 0.02 | 0.03 | 0.03 | 0.05 | 0.06 | 0.03 ± 0.01 |
| Total Quality Score | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 ± 0.00 |

---

### 3.2 Full `lyso1` Dataset (123 Frames: Sweeps 094123, 094348, 100458)

#### Diffracting Frames ($N = 48$)
| Metric | Min | 25th | Median | 75th | 95th | Max | Mean ± Std |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Total Spots ($N_{\text{spots}}$) | 19.00 | 191.50 | 351.50 | 609.25 | 729.10 | 827.00 | 387.27 ± 224.44 |
| Bragg Spots ($N_B$) | 6.00 | 171.75 | 324.00 | 581.75 | 690.20 | 815.00 | 361.71 ± 217.31 |
| Average SNR ($\text{SNR}$) | 20.15 | 133.95 | 219.79 | 452.40 | 1,055.48 | 2,266.78 | 361.77 ± 390.72 |
| Bragg Purity ($P_B$, %) | 31.58 | 88.22 | 92.42 | 96.56 | 98.53 | 98.55 | 91.11 ± 10.02 |
| Resolution Limit ($d_{98}$, Å) | 1.48 | 1.55 | 1.91 | 2.04 | 2.24 | 2.41 | 1.82 ± 0.27 |
| Ice Profile Peak ($I_{\text{ice}}$) | 0.02 | 0.03 | 0.05 | 0.27 | 1.44 | 2.66 | 0.30 ± 0.53 |
| Detected Ice Rings ($N_{\text{ice}}$) | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 ± 0.00 |
| Ice Penalty ($P_{\text{ice}}$) | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.66 | 0.01 ± 0.09 |
| Total Quality Score | 44.56 | 99.58 | 99.86 | 99.89 | 99.92 | 99.94 | 97.87 ± 8.12 |

#### Blank Frames ($N = 75$)
| Metric | Min | 25th | Median | 75th | 95th | Max | Mean ± Std |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Total Spots ($N_{\text{spots}}$) | 1.00 | 2.00 | 2.00 | 2.00 | 2.00 | 3.00 | 2.03 ± 0.23 |
| Bragg Spots ($N_B$) | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 ± 0.00 |
| Average SNR ($\text{SNR}$) | 44,506.38 | 66,817.45 | 66,923.11 | 66,995.88 | 67,054.32 | 140,984.84 | 67,090.23 ± 9666.92 |
| Total Quality Score | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 ± 0.00 |

---

### 3.3 `insul-2` Dataset ($N = 8$)

A high-resolution bovine insulin cryo-dataset exhibiting strong Bragg reflections extending to $d_{98} \approx 1.09$ Å, coupled with significant ice ring powder diffraction on several frames.

| Metric | Min | 25th | Median | 75th | 95th | Max | Mean ± Std |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Total Spots ($N_{\text{spots}}$) | 271.00 | 283.00 | 287.50 | 297.75 | 331.55 | 334.00 | 294.88 ± 21.30 |
| Bragg Spots ($N_B$) | 207.00 | 220.50 | 223.50 | 233.75 | 254.90 | 257.00 | 228.50 ± 15.93 |
| Average SNR ($\text{SNR}$) | 72.12 | 81.97 | 96.51 | 104.47 | 131.26 | 143.10 | 97.74 ± 20.81 |
| Bragg Purity ($P_B$, %) | 72.13 | 76.64 | 77.60 | 79.02 | 81.35 | 81.92 | 77.55 ± 2.81 |
| Resolution Limit ($d_{98}$, Å) | 1.09 | 1.09 | 1.10 | 1.12 | 1.13 | 1.13 | 1.11 ± 0.02 |
| Ice Profile Peak ($I_{\text{ice}}$) | 0.70 | 2.63 | 4.40 | 7.78 | 38.14 | 52.21 | 10.57 ± 16.08 |
| Detected Ice Rings ($N_{\text{ice}}$) | 0.00 | 0.00 | 0.00 | 1.00 | 1.00 | 1.00 | 0.38 ± 0.48 |
| Ice Penalty ($P_{\text{ice}}$) | 0.00 | 0.77 | 2.40 | 7.28 | 10.00 | 10.00 | 4.02 ± 3.93 |
| Total Quality Score | 88.86 | 91.72 | 96.91 | 98.00 | 99.02 | 99.08 | 95.00 ± 3.96 |

> [!NOTE]
> On `insul-2`, ice penalties actively modulate the composite score down from pre-penalty raw scores of 98.73–99.40 to 88.86–99.08. This verifies that the post-logistic ice penalty ($P_{\text{ice}} = \text{clamp}(2 N_{\text{ice}} + \max(0, I_{\text{ice}}-2), 0, 10)$) behaves boundedly as specified in ADR 0020.

---

### 3.4 `thaumatin` Dataset ($N = 50$, Sampled Across 500 Frames)

A consistent, well-centered rotation dataset with sharp, moderate-count spots and clean resolution extending beyond 1.25 Å.

| Metric | Min | 25th | Median | 75th | 95th | Max | Mean ± Std |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Total Spots ($N_{\text{spots}}$) | 106.00 | 143.25 | 157.00 | 165.00 | 180.65 | 185.00 | 154.06 ± 17.30 |
| Bragg Spots ($N_B$) | 82.00 | 122.00 | 136.50 | 145.00 | 160.75 | 165.00 | 133.90 ± 17.55 |
| Average SNR ($\text{SNR}$) | 37.05 | 45.13 | 49.45 | 55.70 | 68.43 | 79.86 | 51.48 ± 9.08 |
| Bragg Purity ($P_B$, %) | 75.18 | 84.57 | 87.31 | 89.24 | 92.67 | 94.07 | 86.79 ± 4.11 |
| Resolution Limit ($d_{98}$, Å) | 1.18 | 1.20 | 1.21 | 1.22 | 1.24 | 1.25 | 1.21 ± 0.02 |
| Ice Profile Peak ($I_{\text{ice}}$) | 0.01 | 0.07 | 0.16 | 0.33 | 1.07 | 2.86 | 0.32 ± 0.48 |
| Detected Ice Rings ($N_{\text{ice}}$) | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 ± 0.00 |
| Ice Penalty ($P_{\text{ice}}$) | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.86 | 0.02 ± 0.12 |
| Total Quality Score | 95.00 | 97.00 | 97.32 | 97.78 | 98.18 | 98.60 | 97.31 ± 0.63 |

---

### 3.5 Cross-Dataset Aggregate: All Diffracting Frames ($N = 106$)

Aggregating all genuine crystalline diffraction frames across all three systems:

| Metric | Min | 25th | Median | 75th | 95th | Max | Mean ± Std |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Total Spots ($N_{\text{spots}}$)** | 19.00 | 153.00 | 170.50 | 336.25 | 676.50 | 827.00 | 270.29 ± 188.68 |
| **Bragg Spots ($N_B$)** | 6.00 | 133.25 | 150.50 | 306.75 | 627.25 | 815.00 | 244.20 ± 183.19 |
| **Average SNR ($\text{SNR}$)** | 20.15 | 48.19 | 64.52 | 195.82 | 768.08 | 2,266.78 | 195.48 ± 303.69 |
| **Bragg Purity ($P_B$, %)** | 31.58 | 85.16 | 88.64 | 92.56 | 98.15 | 98.55 | 88.05 ± 8.21 |
| **Resolution Limit ($d_{98}$, Å)** | 1.09 | 1.20 | 1.23 | 1.64 | 2.10 | 2.41 | 1.48 ± 0.36 |
| **Ice Profile Peak ($I_{\text{ice}}$)** | 0.01 | 0.04 | 0.13 | 0.44 | 2.99 | 52.21 | 1.08 ± 5.20 |
| **Detected Ice Rings ($N_{\text{ice}}$)** | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 1.00 | 0.03 ± 0.17 |
| **Ice Penalty ($P_{\text{ice}}$)** | 0.00 | 0.00 | 0.00 | 0.00 | 0.99 | 10.00 | 0.32 ± 1.52 |
| **Total Quality Score** | 44.56 | 97.13 | 97.88 | 99.84 | 99.91 | 99.94 | 97.39 ± 5.64 |

---

## 4. Logit Factor Decomposition and Saturation Analysis

The current scoring formula (implemented in `csrc/mxspots.c` per ADR 0020) defines the logit index $z$ as:
$$z = w_0 + T_N + T_S + T_P + T_{\text{res}}$$
where:
- Baseline Intercept: $w_0 = -6.50$
- Bragg Spot Count Term: $T_N = 0.75 \ln(1 + N_B)$
- Bragg SNR Term: $T_S = 1.20 \ln(1 + \text{SNR})$
- Bragg Purity Term: $T_P = 1.40 \left(\frac{P_B}{100}\right)$
- Resolution Term: $T_{\text{res}} = 0.50 s_{\text{res}}$, where $s_{\text{res}} = \text{clamp}\left(\frac{4.0 - d_{98}}{4.0 - 1.2}, 0.0, 1.0\right)$
- Raw Logistic Score: $\text{Score}_{\text{raw}} = \frac{100}{1 + \exp(-z)}$

### 4.1 Logit Component Distributions Across Diffracting Frames ($N = 106$)

| Component Term | Min | 25th | Median | 75th | 95th | Max | Mean ± Std | Dynamic Spread ($\Delta$) | Share of Spread |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Intercept ($w_0$) | -6.50 | -6.50 | -6.50 | -6.50 | -6.50 | -6.50 | -6.50 ± 0.00 | 0.00 | 0.0% |
| Bragg Count ($T_N$) | 1.46 | 3.67 | 3.77 | 4.30 | 4.83 | 5.03 | 3.95 ± 0.52 | 3.57 | 34.5% |
| Bragg SNR ($T_S$) | 3.66 | 4.67 | 5.02 | 6.34 | 7.97 | 9.27 | 5.59 ± 1.20 | 5.61 | 54.3% |
| Bragg Purity ($T_P$) | 0.44 | 1.19 | 1.24 | 1.30 | 1.37 | 1.38 | 1.23 ± 0.11 | 0.94 | 9.1% |
| Resolution Limit ($T_{\text{res}}$) | 0.28 | 0.42 | 0.49 | 0.50 | 0.50 | 0.50 | 0.45 ± 0.06 | 0.22 | 2.1% |
| **Combined Logit ($z$)** | **-0.22** | **3.56** | **3.89** | **6.47** | **7.04** | **7.50** | **4.72 ± 1.55** | **7.72** | **100.0%** |
| **Raw Score ($\text{Score}_{\text{raw}}$)** | **44.56** | **97.24** | **97.99** | **99.85** | **99.91** | **99.94** | **97.71 ± 5.50** | **55.38** | — |

### 4.2 Why Current Scores Saturated at 97.0–99.9

In logistic regression, the transformation $\sigma(z) = \frac{1}{1 + e^{-z}}$ exhibits sensitivity (high first derivative $\sigma'(z) > 0.10$) only in the linear transition zone $z \in [-2.0, +2.0]$:
- At $z = 0$: $\sigma(0) = 50.0$
- At $z = 1.0$: $\sigma(1) = 73.1$
- At $z = 2.0$: $\sigma(2) = 88.1$
- At $z = 3.0$: $\sigma(3) = 95.3$
- At $z = 4.0$: $\sigma(4) = 98.2$
- At $z = 6.0$: $\sigma(6) = 99.75$

In the empirical data:
1. **Under-calibrated Intercept**: The baseline intercept $w_0 = -6.50$ is overwhelmed by the positive terms. Even for modest diffraction ($N_B = 100, \text{SNR} = 30$), $T_N = 0.75 \ln(101) = 3.46$ and $T_S = 1.20 \ln(31) = 4.12$. Their sum $T_N + T_S = 7.58$ already exceeds $|w_0|$.
2. **Elevated Empirical SNR**: Experimental Pilatus/Eiger detector data routinely yield average Bragg SNRs ranging from 40 to 800+. Consequently, $T_S = 1.20 \ln(1+\text{SNR})$ regularly evaluates to $5.0\text{--}8.0$.
3. **Severe Sigmoidal Plateauing**: With median $z = 3.89$ (and median $z = 6.54$ in `lyso1`), **over 90% of diffracting frames fall into the asymptotic plateau $z > 3.5$**, where $e^{-z} \to 0$. In this regime, large physical differences in crystal quality (e.g., $N_B = 150, \text{SNR} = 45$ vs $N_B = 800, \text{SNR} = 800$) produce identical scores: $97.3$ vs $99.9$.

---

## 5. Root Cause Diagnosis: Blank Frame $\text{SNR} \ge 66,000$ Artifact

### 5.1 Artifact Manifestation

Across blank frames in `lyso1` (e.g. `lyso_1-0616-094123_0001.cbf` through `0007.cbf`), the spot engine reported:
```text
lyso_1-0616-094123_0001.cbf: spot_count=2, bragg_spots=0, avg_snr=66667.16, d_min=999.0, score=0.0
lyso_1-0616-094123_0002.cbf: spot_count=2, bragg_spots=0, avg_snr=66670.70, d_min=999.0, score=0.0
lyso_1-0616-094123_0007.cbf: spot_count=2, bragg_spots=0, avg_snr=66690.79, d_min=999.0, score=0.0
```

### 5.2 Step-by-Step Mechanical Root Cause

1. **Persistent Detector Hot Pixel Cluster**:
   Direct inspection of pixel data across the `lyso1` frames reveals an unmasked defective pixel pair:
   - Pixel $(y=1135, x=2219)$: intensity $\approx 73,800$ ADU
   - Pixel $(y=1135, x=2220)$: intensity $\approx 133,327$ ADU
   This defect persists across all 123 frames regardless of crystal presence.

2. **Annular Background Clamping on Blank Detectors**:
   In `csrc/mxspots.c`, local background is estimated using an annular query around candidate pixels (`sat_query_annulus`, $R_{\text{in}} = 3, R_{\text{out}} = 6$):
   ```c
   double mean = s_ann / n_ann;
   double var = (sq_ann / n_ann) - (mean * mean);
   double std = (var > 0.0) ? sqrt(var) : 1.0;
   *out_mean = (float)mean;
   *out_std = (std < 1.0) ? 1.0f : (float)std;
   ```
   On a blank detector frame, all surrounding pixels in the annulus are zero. Thus, $\text{mean} = 0.0$ and $\text{var} = 0.0$.
   The code clamps $\text{std}$ to the fallback floor of `1.0f`.
   The signal-to-noise ratio for the hot pixel evaluates to:
   $$\text{SNR} = \frac{v - \text{bg}}{\text{std}} = \frac{133,327.0 - 0.0}{1.0} = 133,327.0$$

3. **Regularity Filtering on Diffracting vs Blank Frames**:
   - **On Diffracting Frames**: `mxspots_analyze_regularity()` constructs difference-vector recurrence graphs. Only spots belonging to periodic lattice graphs of component size $\ge 5$ are marked as Bragg reflections. The isolated hot pixel is rejected. Bragg SNR is computed strictly as:
     $$\text{avg\_snr} = \frac{\sum_{i \in \text{Bragg}} \text{SNR}_i}{N_B}$$
     which yields normal values ($20\text{--}800$).
   - **On Blank Frames**: No Bragg spots exist ($N_B = 0$). Only 2 raw spots exist on the entire frame: the hot pixel spot ($\text{SNR} = 133,326.8$) and 1 random noise pixel ($\text{SNR} = 7.5$).

4. **The C Engine Fallback Bug (`csrc/mxspots.c` lines 1028–1034)**:
   In `mxspots_score()`:
   ```c
   out_score->spot_count = spot_count;
   if (out_score->avg_snr <= 0.0f && spots != NULL && spot_count > 0) {
       double sum_snr = 0.0;
       for (int i = 0; i < spot_count; ++i) {
           sum_snr += spots[i].snr;
       }
       out_score->avg_snr = (float)(sum_snr / spot_count);
   }
   ```
   Because $N_B = 0$, `out_score->avg_snr` was initially $0.0\text{f}$.
   The condition `out_score->avg_snr <= 0.0f && spot_count > 0` triggered, causing the engine to average the SNR across **all** raw spots on the frame:
   $$\text{avg\_snr} = \frac{133,326.8 + 7.5}{2} = 66,667.16$$

5. **Structural Leakage**:
   At line 1040, the engine enforces the strict Bragg gate:
   ```c
   if (n_bragg <= 0) {
       out_score->d_min = 999.0f;
       out_score->score = 0.0f;
       return 0;
   }
   ```
   While this successfully resets `out_score->score = 0.0f`, **`out_score->avg_snr` was already modified and is never cleared**. It is returned to Python in `ScoreResult.avg_snr`, corrupting logs and any downstream model that evaluates logit formulas directly.

---

## 6. Recommendations for Issue #91 (Scoring Model Prototyping)

1. **Fix `avg_snr` Zeroing in `csrc/mxspots.c`**:
   Whenever $N_B = 0$, `out_score->avg_snr` must remain strictly `0.0f` (or NaN). The fallback branch that averages non-Bragg candidate spots should be removed or conditioned on $N_B > 0$.

2. **Recalibrate Scoring Function Dynamic Range**:
   To eliminate the 97–99.9 saturation plateau, candidate models in #91 must shift the typical operating range of $z$ from $[3.5, 7.5]$ back down to $[-1.5, +1.5]$:
   - Increase baseline intercept $|w_0|$ (e.g. from $-6.50$ to $-10.0$ or $-12.0$).
   - Dampen the SNR logit multiplier (e.g. from $1.20$ to $0.40\text{--}0.60$) or use a sublinear / soft-saturating transform $\ln(1 + \ln(1 + \text{SNR}))$.
   - Re-anchor the sigmoid mid-point so that an average diffracting frame ($N_B \approx 150, \text{SNR} \approx 65, d_{98} \approx 1.5$ Å) maps to a score of $\approx 70.0$, weak edge diffraction maps to $30\text{--}50$, and exceptional diffraction maps to $90\text{--}98$.
