# 0008 - Difference Vector Recurrence Graph for Indexing-Free Bragg Regularity

To evaluate crystal diffraction quality and detect split crystals without full 3D unit cell indexing, `mxspots` clusters pairwise reciprocal difference vectors and constructs a union-find adjacency graph. This quantifies Bragg regularity (`percentage_regular`) and identifies distinct crystal lattices (`num_lattices`) robustly on sparse, incomplete, or multi-crystal diffraction frames.
