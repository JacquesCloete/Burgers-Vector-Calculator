% Burger_Vector_Plotter_2D

%% Introduction
% Computes the Burgers vector at each point in 2D space from elastic
% strain and lattice rotation data, producing a field of quiver plots.
% Optionally, the program can also plot where it discerns dislocation
% centres to be, which can help in determining dislocation lines from the
% strains array only.

% Data has been provided for an example test material composed of tungsten,
% though data from other experiments can be used as well. However, it falls
% onto the responsibility of the user to properly implement the new data as
% demonstrated here.

% This program can be used in conjunction with Burgers_Vector_Calculator to
% detect and locate dislocations in a material specimen which can then be
% more closely examined using the Calculator.


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

% Ensure that the coordinate system used by your strain tensor is such that
% the x- and y-axes are in the same plane as the 2D field of data!


% Load the data:
% The function below loads the example data; if you want to use your own
% data then instead implement the method of loading it here:
disp('Loading Experimental Data...')
[strains, Xgrid, Ygrid, GND] = experimental_strains_example_init_2D;
disp('Experimental Data Loaded')

% Convert the strains array into the displacement gradient field, which we
% shall use to perform the numerical integration:
beta = strains2beta_2D(strains);

% Define the interval between adjacent voxels of input data in metres:
% The program assumes that the voxels are equally spaced in each direction.
% The value for the example data is 3.9 nm, change to suit your own data 
% accordingly.
interval = 3.9e-9;    % (m)


%% Computing The Field of Burgers Vectors
hold on

% Determine the limits for the x- and y-coordinates across which to plot:
% (Note the upper limit is the first entry for the x-direction)
% Make sure that the first or last values in m or n, plus the extra
% distance equal to the half-width, does not exceed the limits of the data
% array!
m = 4*3.9e-9:interval:48*3.9e-9;    % x-limits (m)
n = 35*3.9e-9:-interval:4*3.9e-9;    % y-limits (m)

% Set up 3D arrays of zeros to house the computed Burgers vectors and their
% directions:
b_computed = zeros(length(n),length(m),3);
b_directions = zeros(length(n),length(m),3);


% Establish the size of the square about which each Burgers circuit will be
% defined. Given a pixel at the centre of the square, half_width defines
% the closest distance in pixels of each of the square's surfaces to the
% centre pixel. For example, if half_width = 2, the square is 5x5 pixels.

% It should be noted that there is a trade-off here; reducing half_width
% will make your dislocation 'tubes' thinner, as fewer Burgers circuits will
% contain the true dislocation lines, but error will be increased as the
% number of points used for the numerical integration is reduced. Accuracy
% is very high until around half_width = 4, and I do not recommend using
% half_width = 1. half_width = 2 has been used to achieve thin dislocation
% tubes at the cost of some accuracy.
half_width = 3;     % (pixels)


% Establish a cutoff value for the calculated Burgers vector magnitude;
% any calculated Burgers vectors with magnitudes below the cutoff are set
% to zero. Increase this to filter out more quivers.
cutoff = 1.5e-10;   % (m)


% Iterate through the chosen limits, computing the Burgers vector at each
% position:
disp('Starting Burgers Circuit Iterations...')
for jloop = 1:numel(n)
    j = n(jloop);
    for iloop = 1:numel(m)
        i = m(iloop);
        % Define the x-limits of the input data in lab coordinates (in m)
        xlims = [i-(interval*half_width) i+(interval*half_width)];
        
        % Define the y-limits of the input data in lab coordinates (in m)
        % (Note the upper limit is the first entry for the y-direction)
        ylims = [j+(interval*half_width) j-(interval*half_width)];
        
        
        % Snippet of relevant data for this iteration:
        beta_snippet = snippet_2D(beta,xlims,ylims,Xgrid(1,1),Ygrid(1,1),interval);
        
        
        % Compute the Burgers vector from this snippet:
        b_computed(jloop,iloop,:) = find_burgers_vector_fast_2D(interval,xlims,ylims,beta_snippet);
        
        
        % If the calculated Burgers vector magnitude is below the
        % cutoff, it is likely just noise and thus set to 0
        % I also had to add an upper bound for this data specifically to
        % ignore the results from the grain boundary
        if norm(squeeze(b_computed(jloop,iloop,:))) <= cutoff || norm(squeeze(b_computed(jloop,iloop,:))) >= 2.5e-10
            b_computed(jloop,iloop,:) = 0;
        end
        
        % Compute the direction of the Burgers vector from this snippet:
        b_directions(jloop,iloop,:) = b_computed(jloop,iloop,:)/norm(squeeze(b_computed(jloop,iloop,:)));
    end
end
disp('Iterations Complete')

% I considered incrementing through each index in the array in one large 
% for-loop rather than across the columns and rows in an effort to hasten
% the iteration process, but it actually turned out to take marginally
% longer overall. If you can think of methods to speed up the process, 
% please feel free to implement them and make mention of them to me.

%% Plotting the Computed Burgers Vectors

% Plot quivers depicting the coordinate axes. Change the position and size
% to suit your own material specimen.
quiver3(-1e-9,-1*3.9e-9,0,50*3.9e-9./0.9,0,0, 'r', 'LineWidth', 1.5);    % x-axis
quiver3(-1e-9,-1*3.9e-9,0,0,50*3.9e-9./0.9,0, 'g', 'LineWidth', 1.5);    % y-axis
quiver3(-1e-9,-1*3.9e-9,0,0,0,1*3.9e-9./0.9, 'b', 'LineWidth', 1.5);    % z-axis


% Produce meshgrids across which to plot the set of computed Burgers vectors
[mmesh,nmesh] = meshgrid(m,n);

% Plot the computed Burgers vectors. These have been scaled up by an
% additional 10 for improved visibility.
quiver3(mmesh,nmesh,zeros(numel(n),numel(m)),b_computed(1:numel(n),1:numel(m),1)*10,b_computed(1:numel(n),1:numel(m),2)*10,b_computed(1:numel(n),1:numel(m),3)*10, 'k');

% Find the reference Burgers vector direction:
b_ref = Reference_Burgers_Vectors_2D;

% Plot the reference Burgers vector direction:
quiver3(6e-9,6e-9,0,b_ref(1)*1e-8,b_ref(2)*1e-8,b_ref(3)*1e-8,'m');

axis equal
view([0 0 1])
camup([0 -1 0])