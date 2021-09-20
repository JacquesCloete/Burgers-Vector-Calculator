function [strains, Xgrid, Ygrid, GND] = experimental_strains_example_init_2D
% experimental_strains_example_init
% Loads an example set of elastic strain and lattice rotation data (2D)

% Data was provided by Hongbing Yu, reference to the associated paper is as
% follows:
% H. Yu et al. “Mapping the full lattice strain tensor of a single dislocation by high angular resolution transmission Kikuchi diffraction (HR-TKD)”. In: Scripta Materialia 164 (2019),pp. 36–41. doi:https://doi.org/10.1016/j.scriptamat.2018.12.039.


load('tkdwithjunliang3_3_data.mat');   % loads the raw data provided

% Extract the X-Grid and Y-Grid and Establish the empty 4D strains array:
Xgrid = Maps.X*3.9e-9;
Ygrid = flip(Maps.Y,1)*3.9e-9;
dimensions = size(Xgrid);
strains = zeros(dimensions(1),dimensions(2),3,3);


% Insert the data into the strains array:
% (Note that the arrays are flipped so that the first entry corresponds to 
% the most positive y-coordinate and least positive x-coordinate)
% NOTE: The strain tensor provided is deviatoric, but the full strain
% tensor is required. The assumption of out-of-plane stress being zero is
% necessary to determine the full strain tensor.
v = 0.28;   % Poisson's Ratio of Tungsten
hydrostatic_strain = -v/(1+v)*(flip(flip(Maps.E11_F,1),2) + flip(flip(Maps.E22_F,1),2) + (1-v)/v*flip(flip(Maps.E33_F,1),2));   % Derived from Hooke's Law

strains(:,:,1,1) = flip(flip(Maps.E11_F,1),2) + hydrostatic_strain;
strains(:,:,1,2) = flip(flip(Maps.E12_F,1),2);
strains(:,:,1,3) = flip(flip(Maps.E13_F,1),2);
strains(:,:,2,2) = flip(flip(Maps.E22_F,1),2) + hydrostatic_strain;
strains(:,:,2,3) = flip(flip(Maps.E23_F,1),2);
strains(:,:,3,3) = flip(flip(Maps.E33_F,1),2) + hydrostatic_strain;

strains(:,:,2,1) = flip(flip(Maps.W21_F1,1),2);
strains(:,:,3,1) = flip(flip(Maps.W13_F1,1),2);
strains(:,:,3,2) = flip(flip(Maps.W32_F1,1),2);