# 0012 - Synthetic Diffraction Frames for Testing

To support rapid, deterministic test-driven development across spot finding, ice masking, and multi-lattice detection without committing massive binary detector geometries (e.g., CBF, HDF5, MCCD) to the repository, `mxspots` procedurally generates 2D synthetic frame arrays on the fly. These are rendered by depositing 2D Gaussian spot profiles and intensity powder rings directly onto an empty frame using coordinates parsed from XDS reference outputs (`SPOT.XDS`).
