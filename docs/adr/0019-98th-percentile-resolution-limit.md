# 0019 - 98th Percentile Resolution Limit

## Context
ADR 0011 established computing frame resolution limits over confirmed Bragg reflections using a 95th percentile metric ($d_{95}$) to eliminate isolated false positives and noise peaks. However, discarding 5% of confirmed Bragg reflections proved overly conservative for high-resolution frames with sparse but genuine reflections near detector edges.

## Decision
Change the resolution limit metric from the 95th percentile ($d_{95}$) to the 98th percentile ($d_{98}$) computed exclusively over confirmed Bragg reflections:
1. Sort confirmed Bragg reflections ascending by $d$-spacing.
2. Select the cutoff index $k = \lfloor 0.02 \times (N_{\text{bragg}} - 1) \rfloor$ so that 98% of reflections satisfy $d \ge d_{98}$.
3. Feed $d_{98}$ into the composite scoring formula via $s_{\text{res}} = \text{clamp}\left(\frac{4.0 - d_{98}}{4.0 - 1.2}, 0.0, 1.0\right)$.
4. Format CLI output as `X.XX Å (98th percentile)`.
