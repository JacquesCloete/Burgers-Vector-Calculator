# Burgers-Vector-Calculator
Program that computes (and plots) Burgers vectors of dislocations entirely from elastic strain and lattice rotation data. Designed and created by Jacques Cloete under supervision from Felix Hofmann and Edmund Tarleton. This is the supplementary program to the following paper:

**Cloete J**, Tarleton E, Hofmann F. (2022) Computation of Burgers Vectors from Elastic Strain and Lattice Rotation Data. *Proc. R. Soc. A*, 478(2263)
DOI: [https://doi.org/10.1098/rspa.2021.0909](https://doi.org/10.1098/rspa.2021.0909)

It is encouraged that the user reads the paper prior to operating the program.

The program has both a 3D and 2D variant, each of which comprised of a Plotter and a Calculator:

The Plotter finds an estimate of the Burgers vector at each point in space and plots them as a field of 3D quivers.

The Calculator allows for manual targeting of a dislocation or group of dislocations to accurately compute the corresponding Burgers vector.

Thus, dislocations identified and located using the Plotter can then be more closely examined using the Calculator.

Both the 3D and 2D versions are configured to use the experimental data sets provided, though the program can be adjusted to take any elastic strain and lattice rotation data set (this implementation is the responsibility of the user, however, and it is suggested that the user observes how the existing data sets were implemented in order to do this).

The following images were produced using the 3D Plotter, loaded with the experimental data set provided and with the microcrystal shape and known dislocation lines superposed onto the plot. The black regions are composed of quivers, i.e. where the Plotter has computed non-negligible Burgers vectors.

<img src="https://github.com/JacquesCloete/Burgers-Vector-Calculator/blob/main/Burgers_Vector_Program_Cloete_J_FINAL/3D/goodfigure1.png" width="221" height="300" align="left" /> <img src="https://github.com/JacquesCloete/Burgers-Vector-Calculator/blob/main/Burgers_Vector_Program_Cloete_J_FINAL/3D/goodfigure2.png" width="284" height="300" align="center" /> <img src="https://github.com/JacquesCloete/Burgers-Vector-Calculator/blob/main/Burgers_Vector_Program_Cloete_J_FINAL/3D/goodfigure3.png" width="264" height="300" align="right" />

For further information on both the 3D and 2D versions of the program, refer to the readme.txt files in their respective folders.


## Zenodo Archive

A Zenodo archive, preserved at the time of publication, is also available with DOI below:

[![DOI](https://zenodo.org/badge/408582523.svg)](https://zenodo.org/badge/latestdoi/408582523)

## Original Datasets Used

BCDI Data:  Hofmann F, Phillips NW, Das S, Karamched P, Hughes GM, Douglas JO, Cha W, Liu W.
2020 Nanoscale imaging of the full strain tensor of specific dislocations extracted from a bulk
sample. Phys. Rev. Materials 4, 013801. <br />
DOI: https://doi.org/10.1103/PhysRevMaterials.4.013801

HR-TKD Data: Yu H, Liu J, Karamched P, Wilkinson AJ, Hofmann F. 2019 Mapping the full lattice strain tensor
of a single dislocation by high angular resolution transmission Kikuchi diffraction (HR-TKD).
Scripta Materialia 164, 36–41. <br />
DOI: https://doi.org/10.1016/j.scriptamat.2018.12.039

## License

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
