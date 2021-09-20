% This script was made to find the standard deviation for the noise of the 
% 2D data set to be used for analysis in the program's associated paper

% Data was provided by Hongbing Yu, reference to the associated paper is as
% follows:
% H. Yu et al. “Mapping the full lattice strain tensor of a single dislocation by high angular resolution transmission Kikuchi diffraction (HR-TKD)”. In: Scripta Materialia 164 (2019),pp. 36–41. doi:https://doi.org/10.1016/j.scriptamat.2018.12.039.

load('tkdwithjunliang3_3_data.mat');


E11 = Maps.E11_F(12:26,34:51);
E12 = Maps.E12_F(19:38,34:51);
E13 = Maps.E13_F(26:38,26:51);
E22 = Maps.E22_F(1:20,31:51);
E23 = Maps.E23_F(1:20,31:51);
E33 = Maps.E33_F(1:38,34:51);

W21 = Maps.W21_F1(1:20,34:51);
W13 = Maps.W13_F1(26:38,26:51);
W32 = Maps.W32_F1(1:25,34:51);

stdev = zeros(1,9);

stdev(1) = std(E11,0,'all');
stdev(2) = std(E12,0,'all');
stdev(3) = std(E13,0,'all');
stdev(4) = std(E22,0,'all');
stdev(5) = std(E23,0,'all');
stdev(6) = std(E33,0,'all');
stdev(7) = std(W21,0,'all');
stdev(8) = std(W13,0,'all');
stdev(9) = std(W32,0,'all');

stev_mean = mean(stdev,'all');