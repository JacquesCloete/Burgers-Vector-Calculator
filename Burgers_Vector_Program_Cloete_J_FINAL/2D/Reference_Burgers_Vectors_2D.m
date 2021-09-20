function b_ref = Reference_Burgers_Vectors_2D
% Reference_Burgers_Vectors_2D
% Calculates the direction of the Burgers vector given the direction in
% crystal lattice coordinates and the relative orientation between the
% crystal lattice and lab coordinates

% Data was provided by Hongbing Yu, reference to the associated paper is as
% follows:
% H. Yu et al. “Mapping the full lattice strain tensor of a single dislocation by high angular resolution transmission Kikuchi diffraction (HR-TKD)”. In: Scripta Materialia 164 (2019),pp. 36–41. doi:https://doi.org/10.1016/j.scriptamat.2018.12.039.


phi1 = 227.8;
phi2 = 23.2;
phi3 = 322.3;

R = [cosd(phi1)*cosd(phi3)-cosd(phi2)*sind(phi1)*sind(phi3) -cosd(phi1)*sind(phi3)-cosd(phi2)*cosd(phi3)*sind(phi1) sind(phi1)*sind(phi2);
     cosd(phi3)*sind(phi1)+cosd(phi1)*cosd(phi2)*sind(phi3) cosd(phi1)*cosd(phi2)*cosd(phi3)-sind(phi1)*sind(phi3) -cosd(phi1)*sind(phi2);
     sind(phi2)*sind(phi3) cosd(phi3)*sind(phi2) cosd(phi2)];

b_lattice = [1; -1; 1];

b_ref = -R*b_lattice;