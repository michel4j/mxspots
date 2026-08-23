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
The standalone C shared library (`.so`/`.dylib`/`.dll`) dynamically loaded via `ctypes` that performs memory-direct spot finding and indexing.
_Avoid_: C++ runtime, Python extension module

**Spot Finding**:
The process of detecting Bragg spots within a frame array using the Spot Engine according to configured thresholds and background parameters.
_Avoid_: Peak picking, blob detection

**Ice Ring**:
A concentric powder diffraction ring arising from vitreous/hexagonal water ice crystal contamination at characteristic Bragg $d$-spacings (e.g. 3.90 Å, 3.67 Å, 3.44 Å, 2.25 Å, 2.07 Å, 1.92 Å).
_Avoid_: Powder ring, background artifact

**Ice Ring Masking**:
The automatic detection and exclusion of pixels or candidate spots falling within contaminated ice ring resolution shells during spot finding and quality scoring.
_Avoid_: Ring stripping, peak cutting

**Resolution Limit**:
The highest resolution (smallest $d$-spacing in Angstroms, evaluated at the 95th percentile $d_{95}$) at which statistically significant diffraction spots are identified on a frame.
_Avoid_: Max resolution, high-resolution cutoff

**Bragg Regularity**:
The percentage of detected spots on a frame that belong to periodic reciprocal lattice recurrence graphs.
_Avoid_: Spot periodicity, lattice score, regularity rate

**Multi-Lattice / Split Crystal**:
The presence of multiple independent crystal lattices detected on a single diffraction frame, identified by distinct connected components in the difference vector recurrence graph.
_Avoid_: Overlapping frames, double lattice, multiple crystals

**Indexing**:
The assignment of 3D reciprocal lattice vectors to observed spots using the Spot Engine's FFT routine to determine unit cell parameters.
_Avoid_: Auto-indexing, lattice mapping

**Percentage Indexed**:
The percentage of detected spots on a frame that are successfully fitted to the indexing lattice model.
_Avoid_: Indexing rate, indexed spot fraction

**Composite Quality Score**:
A normalized 0–100 quality metric computed by the Spot Engine integrating spot count, Bragg regularity / percentage indexed, average SNR, and 95th percentile resolution limit ($d_{95}$) with penalties for ice contamination and multi-lattice split crystals.
_Avoid_: Frame score, quality index, total score

**XDS Spot File**:
A reference file (e.g. `SPOT.XDS`) containing spot coordinates $(x, y, z)$ and intensities used by test fixtures to synthesize test frames.
_Avoid_: Spot list file

**Synthetic Frame**:
A programmatically generated 2D frame array created by rendering 2D Gaussian spot profiles from an XDS spot file with noise, cached for fast testing.
_Avoid_: Fake frame, mock image

## Workflows & Commands

**mxspots.findspots**:
Command and Python interface for detecting spots on diffraction frames and outputting spot coordinates and intensities.

**mxspots.score**:
Command and Python interface for computing frame quality metrics including spot count, SNR, resolution limit, ice contamination scores, Bragg regularity, and composite quality score.

**mxspots.index**:
Command and Python interface for running lattice indexing on detected spots and calculating percentage indexed.
