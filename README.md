# chimes_CGD

These scripts are designed to operate on the training_data.xyzf file, for our new, unpublished ChIMES-C dataset. A description of this dataset can be found in the second tab of Carbon_data_defs.xlsx.

These scripts should be able to run on UM-ARC and LLNL-LC, but have only been tested on the latter so far. 

The scripts have also been configured for MPI/parallel runs via sbatch.

To run through everything on ARC-TS, execute:

Step 1: ./get_clusters.sh training_data.xyzf

- Expects ~/codes/chimes_lsq-myLLfork/modfiles/UM-ARC.mod to exist
- This depends on helpers.py and extract clusters.cpp
- This launches a series of single-node (36 proc) jobs
- This script produces lists of 2, 3, and 4-body clusters in both "r" and "s" style transformations for every frame in the training trajectory


Step 2: get_histograms.sh

- Expects ~/codes/chimes_lsq-myLLfork/modfiles/UM-ARC.mod to exist
- This depends on calc_cluster_distance_histograms-mpi.cpp
- This launches a series of single-node (36 proc) jobs
- This script produces cluster distance histograms for each framein both "r" and "s" style transformations for every frame in the training trajectory -- referred to as "d" in CluUQ powerpoint

Step 3: get_cluster_similarities.sh

- This depends on CluUQ_similarity.cpp
- Computes a scalar "distance" (called "\Delta_{CG}" in our powerpoint) between every set of two frames
- This script produces dist_matrix.dat, which can be immediately plotted with gnuplot in heatmap mode (splot/pm3d) to visualize frame dissimilarities.
