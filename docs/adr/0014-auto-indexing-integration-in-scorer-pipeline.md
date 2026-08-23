# 0014 - Auto-Indexing Integration in Scorer Pipeline

To provide a comprehensive evaluation of crystal diffraction quality, the `mxspots.score` Python pipeline automatically performs reciprocal space lattice indexing on candidate spots via the Engine's FFT routine. Automatically attempting to index the detected spot list yields an objective `percentage_indexed` rate that heavily weights the final composite quality score, validating whether the detected Bragg regularity corresponds to a physically sound 3D unit cell.
