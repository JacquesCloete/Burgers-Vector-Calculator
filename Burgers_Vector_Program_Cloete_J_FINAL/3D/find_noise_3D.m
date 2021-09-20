% This script was made to find the standard deviation for the noise of the 
% 3D data set to be used for analysis in the program's associated paper

% Data was provided by Felix Hofmann, reference to the associated paper is
% as follows:
% F. Hofmann et al. “Nanoscale imaging of the full strain tensor of specific dislocations extracted from a bulk sample”. In: Phys. Rev. Materials4.1 (2020), p. 013801. doi:https://doi.org/10.1103/PhysRevMaterials.4.013801.


disp('Loading Experimental Data...')
[strains, Xgrid, Ygrid, Zgrid, Sam_red] = experimental_strains_example_init();
disp('Experimental Data Loaded')

interval = 5e-9;

% Define the x-limits of the input data in lab coordinates (in m)
xlims = [-82.5e-9 82.5e-9];    % x-limits (m)

% Define the y-limits of the input data in lab coordinates (in m)
% (Note the upper limit is the first entry for the y-direction)
ylims = [2.5e-9 -152.5e-9];    % y-limits (m)

% Define the z-limits of the input data in lab coordinates (in m)
zlims = [-152.5e-9 2.5e-9];    % z-limits (m)

strains = snippet(strains,xlims,ylims,zlims,Xgrid(1,1,1),Ygrid(end,1,1),Zgrid(1,1,1),interval);

stdev = zeros(1,9);

stdev(1) = std(strains(:,:,:,1,1),0,'all');
stdev(2) = std(strains(:,:,:,1,2),0,'all');
stdev(3) = std(strains(:,:,:,1,3),0,'all');
stdev(4) = std(strains(:,:,:,2,2),0,'all');
stdev(5) = std(strains(:,:,:,2,3),0,'all');
stdev(6) = std(strains(:,:,:,3,3),0,'all');
stdev(7) = std(strains(:,:,:,2,1),0,'all');
stdev(8) = std(strains(:,:,:,3,1),0,'all');
stdev(9) = std(strains(:,:,:,3,2),0,'all');

stev_mean = mean(stdev,'all');