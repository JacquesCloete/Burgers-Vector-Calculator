# Burgers-Vector-Calculator
Program that computes (and plots) Burgers vectors of dislocations entirely from elastic strain and lattice rotation data. Designed and created by Jacques Cloete, with help from Felix Hofmann and Edmund Tarleton. This program is associated with the following paper:

"Computation of Burgers Vectors from Elastic Strain and Lattice Rotation Data"
Jacques Cloete, Felix Hofmann and Edmund Tarleton
(publication pending)

It is encouraged that the user reads the paper prior to operating the program.

The program has both a 3D and 2D variant, each of which comprised of a Plotter and a Calculator:

The Plotter finds an estimate of the Burgers vector at each point in space and plots them as a field of 3D quivers.

The Calculator allows for manual targeting of a dislocation or group of dislocations to accurately compute the corresponding Burgers vector.

Thus, dislocations identified and located using the Plotter can then be more closely examined using the Calculator.

Both the 3D and 2D versions are configured to use the experimental data sets provided, though the program can be adjusted to take any elastic strain and lattice rotation data set (this implementation is the responsibility of the user, however, and it is suggested that the user observes how the existing data sets were implemented in order to do this).

The following images were produced using the 3D Plotter, loaded with the experimental data set provided and with the microcrystal shape and known dislocation lines superposed onto the plot. The black regions are composed of quivers, i.e. where the Plotter has computed non-negligible Burgers vectors.

<img src="https://github.com/JacquesCloete/Burgers-Vector-Calculator/blob/main/Burgers_Vector_Program_Cloete_J_FINAL/3D/goodfigure1.png" width="221" height="300" align="left" /> <img src="https://github.com/JacquesCloete/Burgers-Vector-Calculator/blob/main/Burgers_Vector_Program_Cloete_J_FINAL/3D/goodfigure2.png" width="284" height="300" align="center" /> <img src="https://github.com/JacquesCloete/Burgers-Vector-Calculator/blob/main/Burgers_Vector_Program_Cloete_J_FINAL/3D/goodfigure3.png" width="264" height="300" align="right" />

For further information on both the 3D and 2D versions of the program, refer to the readme.txt files in their respective folders.
