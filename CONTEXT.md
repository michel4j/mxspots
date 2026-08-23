# mxspots Context

`mxspots` provides fast spot analysis and quality assessment for Macromolecular Crystallography (MX) diffraction frames.

## Language

**Diffraction Frame**:
A 2D detector pixel array representing X-ray diffraction intensities from a crystal, loaded via `mxio`.
_Avoid_: Raw image, detector shot

**Spot**:
A localized region of high diffraction intensity on a frame corresponding to a Bragg reflection, defined by `(x, y, d-spacing, intensity)`.
_Avoid_: Peak, reflection, blob

**Spot Engine**:
The standalone C shared library (`.so`/`.dylib`/`.dll`) dynamically loaded via `ctypes` that performs memory-direct spot finding, regularity analysis, and quality scoring.
_Avoid_: C++ runtime, Python extension module

**Spot Finding**:
The process of detecting Bragg spots within a frame array using the Spot Engine according to configured thresholds and background parameters.
_Avoid_: Peak picking, blob detection

**Detector Geometry**:
Metadata describing the physical layout of the detector during data collection, including beam center, distance, wavelength, and pixel size, essential for accurate resolution ($d$-spacing) coordinates and correct SPOT.XDS extraction.
_Avoid_: Detector config, setup parameters

**1D Radial Profile**:
An azimuthally-integrated projection of the 2D frame plotting average pixel intensity against radial distance from the beam center, used for canonical ice ring detection.
_Avoid_: Intensity graph, radial binning

**Ice Ring**:
A concentric powder diffraction ring arising from vitreous/hexagonal water ice crystal contamination at characteristic Bragg $d$-spacings (e.g. 3.90 Å, 3.67 Å, 3.44 Å, 2.25 Å, 2.07 Å, 1.92 Å).
_Avoid_: Powder ring, background artifact

**Ice Ring Masking**:
The automatic detection and exclusion of pixels or candidate spots falling within contaminated ice ring resolution shells during spot finding and quality scoring.
_Avoid_: Ring stripping, peak cutting

**Resolution Limit**:
The highest resolution (smallest $d$-spacing in Angstroms, evaluated at the 95th percentile $d_{95}$) at which statistically significant diffraction spots are identified on a frame.
_Avoid_: Max resolution, high-resolution cutoff

**Bragg Spots**:
The subset of detected spots that belong to regular crystalline periodic lattice graphs (connected components of size $\ge 5$ in reciprocal difference-vector recurrence space).
_Avoid_: Regular spots, lattice points

**Bragg Percent**:
The percentage of total detected spots that are classified as regular Bragg spots ($100 \times N_{\text{bragg}} / N_{\text{spots}}$).
_Avoid_: Regularity rate, percentage regular, indexed fraction

**Average Bragg Intensity**:
The mean integrated intensity across all regular Bragg spots on a frame.
_Avoid_: Mean intensity, total spot intensity

**Multi-Lattice / Split Crystal**:
The presence of multiple independent crystal lattices detected on a single diffraction frame, identified by distinct connected components in the difference vector recurrence graph.
_Avoid_: Overlapping frames, double lattice, multiple crystals

**Composite Quality Score**:
A normalized 0–100 quality metric computed by the Spot Engine using a Hybrid Gated-Logistic model combining Bragg spot count, Bragg percent, Bragg average intensity, average SNR, and 95th percentile resolution limit ($d_{95}$) with ice contamination penalties. A strict zero-gate applies if no Bragg spots are found ($N_{\text{bragg}} = 0$).
_Avoid_: Frame score, quality index, total score

**XDS Spot File**:
A reference file (e.g. `SPOT.XDS`) containing spot coordinates $(x, y, z)$ and intensities used by test fixtures to synthesize test frames.
_Avoid_: Spot list file

**Synthetic Frame**:
A programmatically generated 2D frame array created by rendering 2D Gaussian spot profiles from an XDS spot file with noise, cached for fast testing.
_Avoid_: Fake frame, mock image

## Workflows & Commands

**mxspots.findspots**:
Command and Python interface for detecting spots on diffraction frames and outputting spot coordinates, intensities, and ice ring summaries.

**mxspots.score**:
Command and Python interface for computing frame quality metrics including spot count, Bragg spots, Bragg percent, average intensity, SNR, resolution limit ($d_{95}$), ice contamination score, number of lattices, and composite quality score.
