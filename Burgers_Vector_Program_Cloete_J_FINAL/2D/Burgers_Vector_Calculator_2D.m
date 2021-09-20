% Burger_Vector_Calculator_2D

%% Introduction
% Computes the Burgers vector for a specified Burgers circuit from elastic
% strain and lattice rotation data and outputs its direction and magnitude.
% Optionally, there is implementation for drawing several similar Burgers
% circuits and finding the mean-average Burgers vector (this could be
% useful if the data is too noisy to provide accurate results from a single
% Burgers circuit).

% Data has been provided for an example test material composed of tungsten,
% though data from other experiments can be used as well. However, it falls
% onto the responsibility of the user to properly implement the new data as
% demonstrated here.

% The reference Burgers vectors for the dislocations found in the example
% material specimen are also provided for comparison.

% This program can be used in conjunction with Burgers_Vector_Plotter to
% 'target' and examine dislocations found in the Plotter in greater detail.


%% Loading the Experimental Data
% Note that the data for the elastic strains and lattice rotations must be
% input into a 4D strains array.

% The strains tensor should be assembled as a
% (y_length) x ((x_length) x 3 x 3 array:

%                              [ e11 e12 e13 ]
%  strains(y_pos,x_pos,:,:) =  [ wz  e22 e23 ]
%                              [ wy  wx  e33 ]

% Where each of the above elements have the definitions found in HLT.

% In addition, when constructing the strains array, the first elements in
% the x-direction must correspond to the most negative x-coordinate, but 
% the first element in the y-direction must correspond to the most positive
% y-coordinate. This is so that the strains array is arranged to imitate a 
% traditional 2D Cartesian coordinate system for the purposes of defining
% the Burgers circuit for numerical integration.

% Note that the most negative x-coordinate and the most positive
% y-coordinate for the sample data are also required for the snippet_2D
% function.


% Load the data:
% The function below loads the example data; if you want to use your own
% data then instead implement the method of loading it here:
disp('Loading Experimental Data...')
[strains, Xgrid, Ygrid,GND] = experimental_strains_example_init_2D;
disp('Experimental Data Loaded')

% NOTE: If the data has already been loaded into your workspace in the
% correct format (e.g. by already having run this program), I suggest
% commenting out the above three 3 lines to significantly hasten runtime.

% Convert the strains array into the displacement gradient field, which we
% shall use to perform the numerical integration:
beta = strains2beta_2D(strains);

% Define the interval between adjacent pixels of input data in metres:
% The program assumes that the pixels are equally spaced in each direction.
% The value for the example data is 3.9nm, change to suit your own data 
% accordingly.
interval = 3.9e-9;    % (m)


%% Defining the Burgers Circuit
% Here we define the values of x and y that determine the square around
% which the Burgers circuit is defined. The larger the ranges, and the
% further the dislocation centre from one of the limits, the more accurate
% the resulting Burgers vector. However it should be noted that fluctuation
% in the elements of the computed Burgers vector is small until the
% lengths of the Burgers circuit are less than 9 pixels and/or the
% dislocation centre is within a few pixels of the Burgers circuit.

% If you want to calculate a set of Burgers vectors to find a mean-average,
% note that the additional Burgers circuits will all be smaller than the
% first, 'outermost' one, so make sure that the limits defined below are
% large and/or the value of i is low. Improved results over a single loop
% is only expected if the data is very noisy.

% Define the x-limits of the input data in lab coordinates (in m)
xlims = [22*3.9e-9 38*3.9e-9];    % x-limits (m)

% Define the y-limits of the input data in lab coordinates (in m)
% (Note the upper limit is the first entry for the y-direction)
ylims = [31*3.9e-9 15*3.9e-9];    % y-limits (m)


%% Computing the Burgers Vector
% We shall extract the required snippet of data for the chosen Burgers
% circuit, and then integrate numerically to compute the Burgers vector.

% Snippet of relevant data for this Burgers Circuit:
beta_snippet = snippet_2D(beta,xlims,ylims,Xgrid(1,1),Ygrid(1,1),interval);

% Determine whether a single Burgers circuit or many similar ones is to be
% used:
average = 'N';

if average == 'N'
    % Compute the Burgers vector from this snippet:
    b_computed = find_burgers_vector_accurate_2D(interval,xlims,ylims,beta_snippet);
    
elseif average == 'Y'
    % Alternatively, compute a set of Burgers vectors for many concentric Burgers
    % circuits from this snippet:
    
    % Number of ADDITIONAL concentric Burgers circuits:
    i = 1;
    
    % Compute the set of Burgers vectors from this snippet:
    b_computed_set = find_burgers_vector_repeated_2D(interval,xlims,ylims,beta_snippet,i);
    
    % Remove outliers:
    for row = 1:3
        [B,TF] = rmoutliers(b_computed_set(row,:),'mean','ThresholdFactor',1.5);    % Adjust threshold factor from 1.5 as desired
        b_computed_set = b_computed_set(:,not(TF));
    end
    
    % Compute the mean-average Burgers vector:
    b_computed = sum(b_computed_set,2)/((i+1));
    
end


% For better insight, we can describe the computed Burgers vector in terms
% of its direction and magnitude:

% Find the unit vector in the direction of the Burgers vector:
b_direction = b_computed/norm(b_computed);

% Find the magnitude of the Burgers vector:
b_magnitude = norm(b_computed);

% Optional - Calculate the direction of the reference Burgers vector:
b_ref = Reference_Burgers_Vectors_2D;
b_ref = b_ref/norm(b_ref);