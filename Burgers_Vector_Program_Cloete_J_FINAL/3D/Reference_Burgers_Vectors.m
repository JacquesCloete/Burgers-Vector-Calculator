function [bv1n, bv2n, bv3n, bv4n, bv5n] = Reference_Burgers_Vectors(Sam_red)
% Reference_Burgers_Vectors
% Calculates the Burgers vectors of each of the dislocations from the
% example data according to the reflection orientations in lab coordinates 
% and observed Burgers vectors in crystal lattice vector notation using g.b
% contrast.

% The unit vectors in the directions of the Burgers vectors calculated here
% should be compared to those calculated using Burgers_Vector_Calculator.

% The following code and data was provided by Felix Hofmann, reference to 
% the associated paper is as follows:
% F. Hofmann et al. “Nanoscale imaging of the full strain tensor of specific dislocations extractedfrom a bulk sample”. In: Phys. Rev. Materials4.1 (2020), p. 013801. doi:https://doi.org/10.1103/PhysRevMaterials.4.013801.


a0 = 3.165*10^-10; %tungsten lattice constant in m

% work out dislocation burgers vectors
% compute crystal basis vectors: 
% reflection order:
%10-1,
%-10-1
%1-10
%0-1-1
%110
%01-1
%using the q directions work out the directions of crystal 100, 010 and 001 axes. Call these n100, n010, n001 respectively. 
n100 = Sam_red(1,1).q_sam - Sam_red(2,1).q_sam;  n100  = n100 ./norm(n100);
n100a = Sam_red(3,1).q_sam + Sam_red(5,1).q_sam; n100a = n100a./norm(n100a); %compute as a check...
n010 = Sam_red(5,1).q_sam - Sam_red(3,1).q_sam;  n010  = n010 ./norm(n010);
n010a = Sam_red(6,1).q_sam - Sam_red(4,1).q_sam; n010a = n010a./norm(n010a); %compute as a check...
n001 = -Sam_red(1,1).q_sam - Sam_red(2,1).q_sam;  n001  = n001 ./norm(n001);
n001a = -Sam_red(4,1).q_sam - Sam_red(6,1).q_sam; n001a = n001a./norm(n001a); %compute as a check...
n001b = cross(n100, n010); % just to check this 001 cross 010 gives 001...

% bv1n to bv5n are the unit normals in the directions of the Burgers
% vectors of dislocations 1 to 5 (use these for comparison):
bv1_dir = -0.5 * [1 1 1]';   bv1 = (bv1_dir(1,1)* n100 + bv1_dir(2,1)* n010 + bv1_dir(3,1)* n001)* a0; bv1n = -bv1./norm(bv1);
bv2_dir = 0.5 * [-1 1 1]';  bv2 = (bv2_dir(1,1)* n100 + bv2_dir(2,1)* n010 + bv2_dir(3,1)* n001)* a0; bv2n = -bv2./norm(bv2);
bv3_dir = -[1 0 0]';         bv3 = (bv3_dir(1,1)* n100 + bv3_dir(2,1)* n010 + bv3_dir(3,1)* n001)* a0; bv3n = -bv3./norm(bv3);
bv4_dir = -0.5 * [1 1 1]';   bv4 = (bv4_dir(1,1)* n100 + bv4_dir(2,1)* n010 + bv4_dir(3,1)* n001)* a0; bv4n = -bv4./norm(bv4);
bv5_dir = -0.5 * [1 -1 -1]'; bv5 = (bv5_dir(1,1)* n100 + bv5_dir(2,1)* n010 + bv5_dir(3,1)* n001)* a0; bv5n = -bv5./norm(bv5);

