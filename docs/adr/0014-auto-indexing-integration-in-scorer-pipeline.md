# 0014 - Decoupled Reciprocal Lattice Indexing from Scorer Pipeline

To maximize frame triage throughput and maintain high-speed scoring, `mxspots.score` evaluates crystal diffraction quality using difference-vector Bragg regularity (`percentage_regular`) without running reciprocal lattice FFT indexing by default. Lattice indexing remains fully available via `mxspots.index` / `index_spots()` when explicit unit cell parameters or indexed spot fractions are required.
