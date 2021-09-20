function strains = mixed(xlims,ylims,zlims,interval,b,v,x_translate,y_translate,z_translate,a,psi,theta,phi)
% mixed
% Establishes the 5D strains array for the mathematical model of a mixed 
% dislocation with sense and Burgers vector in any orientation (determined 
% by inputs)

% The derivation for the mathematics of this script is found in the
% program's associated paper.


% Determine the x-, y- and z-coordinates for which data is to be generated:
x_length = int16(abs(xlims(2) - xlims(1))/interval + 1);
x = linspace(xlims(1)-x_translate,xlims(2)-x_translate,x_length);

y_length = int16(abs(ylims(2) - ylims(1))/interval + 1);
y = linspace(ylims(1)-y_translate,ylims(2)-y_translate,y_length)';

z_length = int16(abs(zlims(2) - zlims(1))/interval + 1);
z = zeros(1,1,z_length);
z(1,1,:) = linspace(zlims(1)-z_translate,zlims(2)-z_translate,z_length);


% Compute the overall 3D rotation matrix as defined in the program's
% associated paper:
r = [cosd(phi) -sind(phi) 0; sind(phi) cosd(phi) 0; 0 0 1]*[cosd(theta) 0 sind(theta); 0 1 0; -sind(theta) 0 cosd(theta)]*[cosd(psi) -sind(psi) 0; sind(psi) cosd(psi) 0; 0 0 1];


% Extract the required elements of the rotation matrix and its transpose:
r11 = r(1,1);
r12 = r(1,2);
r13 = r(1,3);
r21 = r(2,1);
r22 = r(2,2);
r23 = r(2,3);
r31 = r(3,1);
r32 = r(3,2);
r33 = r(3,3);


rt11 = r11;
rt12 = r21;
rt13 = r31;
rt21 = r12;
rt22 = r22;
rt23 = r32;
rt31 = r13;
rt32 = r23;
rt33 = r33;


% Constant A, which appears in the mathematical models of the infinite 
% straight dislocations, is as defined in the program's associated paper
A = b*sind(a)/(4*pi*(1-v));


% We must now apply the rotation to the displacement gradient field for the
% mixed dislocation, as described in the program's associated paper:

% x -> (rt11*x + rt12*y + rt13*z)
% y -> (rt21*x + rt22*y + rt23*z)

beta_11 = -A*(rt21*x + rt22*y + rt23*z)./((rt11*x + rt12*y + rt13*z).^2+(rt21*x + rt22*y + rt23*z).^2).^2.*((1-2*v)*(rt21*x + rt22*y + rt23*z).^2 + (3-2*v)*(rt11*x + rt12*y + rt13*z).^2);
beta_12 = A*(rt11*x + rt12*y + rt13*z)./((rt11*x + rt12*y + rt13*z).^2+(rt21*x + rt22*y + rt23*z).^2).^2.*((1-2*v)*(rt21*x + rt22*y + rt23*z).^2 + (3-2*v)*(rt11*x + rt12*y + rt13*z).^2);
beta_21 = -A*(rt11*x + rt12*y + rt13*z)./((rt11*x + rt12*y + rt13*z).^2+(rt21*x + rt22*y + rt23*z).^2).^2.*((1-2*v)*(rt11*x + rt12*y + rt13*z).^2 + (3-2*v)*(rt21*x + rt22*y + rt23*z).^2);
beta_22 = -A*(rt21*x + rt22*y + rt23*z)./((rt11*x + rt12*y + rt13*z).^2+(rt21*x + rt22*y + rt23*z).^2).^2.*((1-2*v)*(rt21*x + rt22*y + rt23*z).^2 + (1+2*v)*(rt11*x + rt12*y + rt13*z).^2);
beta_31 = -(rt21*x + rt22*y + rt23*z)./((rt11*x + rt12*y + rt13*z).^2 + (rt21*x + rt22*y + rt23*z).^2)*b*cosd(a)/(2*pi);
beta_32 = (rt11*x + rt12*y + rt13*z)./((rt11*x + rt12*y + rt13*z).^2 + (rt21*x + rt22*y + rt23*z).^2)*b*cosd(a)/(2*pi);


new_beta_11 = r11*(beta_11*rt11+beta_12*rt21)...
            + r12*(beta_21*rt11+beta_22*rt21)...
            + r13*(beta_31*rt11+beta_32*rt21);

new_beta_12 = r11*(beta_11*rt12+beta_12*rt22)...
            + r12*(beta_21*rt12+beta_22*rt22)...
            + r13*(beta_31*rt12+beta_32*rt22);

new_beta_13 = r11*(beta_11*rt13+beta_12*rt23)...
            + r12*(beta_21*rt13+beta_22*rt23)...
            + r13*(beta_31*rt13+beta_32*rt23);
        
new_beta_21 = r21*(beta_11*rt11+beta_12*rt21)...
            + r22*(beta_21*rt11+beta_22*rt21)...
            + r23*(beta_31*rt11+beta_32*rt21);

new_beta_22 = r21*(beta_11*rt12+beta_12*rt22)...
            + r22*(beta_21*rt12+beta_22*rt22)...
            + r23*(beta_31*rt12+beta_32*rt22);

new_beta_23 = r21*(beta_11*rt13+beta_12*rt23)...
            + r22*(beta_21*rt13+beta_22*rt23)...
            + r23*(beta_31*rt13+beta_32*rt23);

new_beta_31 = r31*(beta_11*rt11+beta_12*rt21)...
            + r32*(beta_21*rt11+beta_22*rt21)...
            + r33*(beta_31*rt11+beta_32*rt21);

new_beta_32 = r31*(beta_11*rt12+beta_12*rt22)...
            + r32*(beta_21*rt12+beta_22*rt22)...
            + r33*(beta_31*rt12+beta_32*rt22);

new_beta_33 = r31*(beta_11*rt13+beta_12*rt23)...
            + r32*(beta_21*rt13+beta_22*rt23)...
            + r33*(beta_31*rt13+beta_32*rt23);


% To match the required input form for the Burgers vector calculator, we
% now form the 5D strains array from the rotated displacement gradient
% field.

% Strains and lattice rotation elements:
exx = new_beta_11;
exy = 1/2*(new_beta_21 + new_beta_12);
exz = 1/2*(new_beta_13 + new_beta_31);
eyy = new_beta_22;
eyz = 1/2*(new_beta_32 + new_beta_23);
ezz = new_beta_33;

wx = 1/2*(new_beta_32 - new_beta_23);
wy = 1/2*(new_beta_13 - new_beta_31);
wz = 1/2*(new_beta_21 - new_beta_12);

strains = zeros(x_length,y_length,z_length,3,3);

% Assemble the strains array:
strains(:,:,:,1,1) = exx;
strains(:,:,:,1,2) = exy;
strains(:,:,:,1,3) = exz;
strains(:,:,:,2,1) = wz;
strains(:,:,:,2,2) = eyy;
strains(:,:,:,2,3) = eyz;
strains(:,:,:,3,1) = wy;
strains(:,:,:,3,2) = wx;
strains(:,:,:,3,3) = ezz;