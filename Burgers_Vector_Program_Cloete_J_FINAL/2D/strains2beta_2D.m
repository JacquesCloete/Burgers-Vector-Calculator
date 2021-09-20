function beta = strains2beta_2D(strains)
% strains2beta_2D
% Converts the strains array into the displacement gradient field, which we
% shall use to perform the numerical integration

% Create an empty 4D array in which to store the displacement gradient:
beta = zeros(size(strains));

% We must first extract the 2D lattice rotation and elastic strain fields
% from the 4D array:
e11 = squeeze(strains(:,:,1,1));

e12 = squeeze(strains(:,:,1,2));

e13 = squeeze(strains(:,:,1,3));

e22 = squeeze(strains(:,:,2,2));

e23 = squeeze(strains(:,:,2,3));

e33 = squeeze(strains(:,:,3,3));


wz = squeeze(strains(:,:,2,1));

wy = squeeze(strains(:,:,3,1));

wx = squeeze(strains(:,:,3,2));


% Now we shall form the displacement gradient elements from these:
beta(:,:,1,1) = e11;

beta(:,:,1,2) = e12 - wz;

beta(:,:,1,3) = e13 + wy;

beta(:,:,2,1) = e12 + wz;

beta(:,:,2,2) = e22;

beta(:,:,2,3) = e23 - wx;

beta(:,:,3,1) = e13 - wy;

beta(:,:,3,2) = e23 + wx;

beta(:,:,3,3) = e33;