# 0018 - Negative Intensities as Masks

## Context

Detector processing pipelines often produce diffraction frames with masked regions, such as dead pixels, module gaps, and beam-stops. These regions are traditionally captured in a secondary boolean mask array that must be passed around and applied at execution time.

`mxspots` needs a clean and fast way to recognize these invalid regions without increasing memory usage or slowing down integral image calculations with secondary bounds checks.

## Decision

We use negative pixel intensities (typically `-1.0f` or `-2.0f`) encoded directly into the signal `float32` array to act as implicit masks representing dead or excluded pixels.

## Consequences

- The C engine implicitly tests `if (v <= 0.0f)` or `if (v >= 0.0f)` in hot loops, branching cleanly and avoiding any secondary array index lookups.
- Background integration image accumulation bypasses these negative pixels smoothly.
- The Python API does not need an overloaded argument signature taking mask arrays `findspots(..., mask=mask_array)`; clients merely construct `data[mask] = -1.0` before passing to C.
- Reduces memory bandwidth inside the C engine loops.
