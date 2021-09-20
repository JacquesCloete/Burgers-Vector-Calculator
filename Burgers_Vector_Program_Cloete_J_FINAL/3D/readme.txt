This repository was created by Jacques Cloete in September 2021, and contains code associated 
with the following paper:

"Computation of Burgers Vectors from Strain and Lattice Rotation Data"
Jacques Cloete, Felix Hofmann, Edmund Tarleton
[pending publication]


The Plotter iterates through a chosen set of data and computes an approximate Burgers vector 
at every voxel/pixel position using a Burgers circuit of standard shape and size. The set of 
Burgers vectors can then be plotted as quivers in 3D space, and the geometry of the material 
specimen itself, known dislocation lines, etc. can be superposed onto the plot if 
implemented by the user.

The Calculator allows for manual input of Burgers circuit limits to target a dislocation or 
group of dislocations and accurately compute the Burgers vector therein.

The two halves of the program can be used effectively in conjunction with each other, with the 
Plotter able to find dislocations and give estimates of their Burgers vectors and the 
Calculator then being used to investigate the discovered dislocations in further detail. The 
Calculator may also be configured to iterate through several concentric Burgers circuits to 
remove outliers and find an average computed Burgers vector.

The program can be adjusted to suit a wide variety of data; all it requires are the elastic 
strain and lattice rotation fields (in the form of evenly-spaced data arrays). The 
implementation of the data will naturally have to be adjusted depending on how the data is 
originally stored.

As well as working with experimental data, the program features a mathematical model for an 
infinite straight mixed dislocation with sense and Burgers vector being able to take any 
direction as described in the associated paper.

The "goodfigure" figure and the "Dislocations" video are outputs from the Plotter, and should 
be achieved by running the program as-is.

The experimental data as well as the useful "strain_rot_tensor" video was provided by Felix 
Hofmann, reference to the associated paper is as follows:
F. Hofmann et al. “Nanoscale imaging of the full strain tensor of specific dislocations extracted from a bulk sample”. In: Phys. Rev. Materials4.1 (2020), p. 013801. doi:https://doi.org/10.1103/PhysRevMaterials.4.013801.


To run the program, see the Calculator and Plotter scripts. They should explain how the 
program works and how to perform the computations on a step-by-step level. The functions etc. 
have also been annotated to explain the procedure clearly.

The "Test" and "find_noise" scripts have been included in case the user wants to 
achieve\verify the findings used for analysis in the associated paper.

This is the 3D variant of the program, see the "2D" folder for the 2D variant.