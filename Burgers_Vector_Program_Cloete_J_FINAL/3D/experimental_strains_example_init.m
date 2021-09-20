function [strains, Xgrid, Ygrid, Zgrid, Sam_red] = experimental_strains_example_init()
% experimental_strains_example_init
% Loads an example set of elastic strain and lattice rotation data
% For the optional plotting of the material shape, the meshgrids of the x-,
% y- and z- coordinates and all of 'Sam_red' are also required

% Data was provided by Felix Hofmann, reference to the associated paper is
% as follows:
% F. Hofmann et al. “Nanoscale imaging of the full strain tensor of specific dislocations extracted from a bulk sample”. In: Phys. Rev. Materials4.1 (2020), p. 013801. doi:https://doi.org/10.1103/PhysRevMaterials.4.013801.


load('recD_sam_str_rot.mat');   % loads the raw data provided

eps_waf = Sam_red.eps_waf;  % loads the elastic strain tensor data
rot_waf = Sam_red.rot_waf;  % loads the lattice rotation tensor data


% Flip the arrays so that the first entry in the y-column corresponds to
% the most positive y-coordinate (as required for the integration)
strains = flip(eps_waf,1);
rot_waf = flip(rot_waf,1);


% Extract the relevant elements from the rot_waf array and reform them into
% the required 5D array format
strains(:,:,:,2,1) = rot_waf(:,:,:,3);
strains(:,:,:,3,1) = rot_waf(:,:,:,2);
strains(:,:,:,3,2) = rot_waf(:,:,:,1);