% Burger_Vector_Plotter

%% Introduction
% Computes the Burgers vector at each point in 3D space from elastic
% strain and lattice rotation data, producing a field of quiver plots.
% Optionally, the actual dislocation lines and the 3D shape of the material
% itself can be superposed onto this for useful visualisation/insight. An
% M-PEG4 video looking around the figure can also be produced.
% The program can also optionally plot where it discerns dislocation
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
% input into a 5D strains array.

% The strains tensor should be assembled as a
% (y_length) x ((x_length) x (z_length) x 3 x 3 array:

%                                    [ e11 e12 e13 ]
%  strains(y_pos,x_pos,z_pos,:,:) =  [ wz  e22 e23 ]
%                                    [ wy  wx  e33 ]

% Where each of the above elements have the definitions found in HLT.

% In addition, when constructing the strains array, the first elements in
% the x- and z-directions must correspond to the most negative x- and
% z-coordinates, but the first element in the y-direction must correspond
% to the most positive y-coordinate. This is so that the strains array
% is arranged to imitate a traditional 3D Cartesian coordinate system for
% the purposes of defining the Burgers circuit for numerical integration.

% Note that the most negative x- and z-coordinates and the most positive
% y-coordinate for the sample data are also required.


% Load the data:
% The function below loads the example data; if you want to use your own
% data then instead implement the method of loading it here:
disp('Loading Experimental Data...')
[strains, Xgrid, Ygrid, Zgrid, Sam_red] = experimental_strains_example_init();
disp('Experimental Data Loaded')

% Convert the strains array into the displacement gradient field, which we
% shall use to perform the numerical integration:
beta = strains2beta(strains);

% Define the interval between adjacent voxels of input data in metres:
% The program assumes that the voxels are equally spaced in each direction.
% The value for the example data is 5e-9 metres, change to suit your own
% data accordingly.
interval = 5e-9;    % (m)


%% Computing The Field of Burgers Vectors
hold on

% Determine the limits for the x-, y- and z-coordinates across which to
% plot:
m = -282.5e-9:interval:282.5e-9;    % x-limits (m)
% (Note the upper limit is the first entry for the y-direction)
n = 402.5e-9:-interval:-392.5e-9;    % y-limits (m)
o = -282.5e-9:interval:282.5e-9;    % z-limits (m)

% Set up 4D arrays of zeros to house the computed Burgers vectors and their
% directions:
b_computed = zeros(length(n),length(m),length(o),3);
b_directions = zeros(length(n),length(m),length(o),3);


% Establish the size of the cube about which each Burgers circuit will be
% defined. Given a voxel at the centre of the cube, half_width defines
% the closest distance in voxels of each of the cube's surfaces to the
% centre voxel. For example, if half_width = 2, the cube is 5x5x5 voxels.

% It should be noted that there is a trade-off here; reducing half_width
% will make your dislocation 'tubes' thinner, as fewer Burgers circuits will
% contain the true dislocation lines, but error will be increased as the
% number of points used for the numerical integration is reduced. Accuracy
% is very high until around half_width = 4, and I do not recommend using
% half_width = 1. half_width = 2 has been used to achieve thin dislocation
% tubes at the cost of some accuracy.
half_width = 2;     % (voxels)


% Establish a cutoff value for the calculated Burgers vector magnitude;
% any calculated Burgers vectors with magnitudes below the cutoff are set
% to zero. Increase this to filter out more quivers.
cutoff = 1.5e-10;   % (m)


% Iterate through the chosen limits, computing the Burgers vector at each
% position:
% (NOTE : the iteration process is timed and the percentage completion is
% also provided)
disp('Starting Burgers Circuit Iterations...')
tic
for jloop = 1:numel(n)
    j = n(jloop);
    for iloop = 1:numel(m)
        i = m(iloop);
        for kloop = 1:numel(o)
            k = o(kloop);
            % The 'if' has been implemented to limit the region of interest
            % to a cylinder, which roughly corresponds to the shape of the
            % example material specimen. For more general cases this will 
            % be irrelevant, so remove as appropriate (or alternatively,
            % create your own restrictions to better match your material).
            if i^2 + k^2 <= (280e-9)^2
                % Define the x-limits of the input data in lab coordinates (in m)
                xlims = [i-(interval*half_width) i+(interval*half_width)];
                
                % Define the y-limits of the input data in lab coordinates (in m)
                % (Note the upper limit is the first entry for the y-direction)
                ylims = [j+(interval*half_width) j-(interval*half_width)];
                
                % Define the z-limits of the input data in lab coordinates (in m)
                zlims = [k-(interval*half_width) k+(interval*half_width)];
                
                
                % Snippet of relevant data for this iteration:
                beta_snippet = snippet(beta,xlims,ylims,zlims,Xgrid(1,1,1),Ygrid(end,1,1),Zgrid(1,1,1),interval);
                
                
                % Compute the Burgers vector from this snippet:
                
                b_computed(jloop,iloop,kloop,:) = find_burgers_vector_fast(interval,xlims,ylims,zlims,beta_snippet);
                
                
                % If the calculated Burgers vector magnitude is below the
                % cutoff, it is likely just noise and thus set to 0
                if norm(squeeze(b_computed(jloop,iloop,kloop,:))) <= cutoff
                    b_computed(jloop,iloop,kloop,:) = 0;
                end
                
                % Compute the direction of the Burgers vector from this snippet:
                b_directions(jloop,iloop,kloop,:) = b_computed(jloop,iloop,kloop,:)/norm(squeeze(b_computed(jloop,iloop,kloop,:)));
            end
        end
    end
    toc
    
    % Establish percentage completion:
    progress = jloop/numel(n)*100;
    fprintf('%.1f%%\n', progress)
end
disp('Iterations Complete')

% I considered incrementing through each index in the array in one large 
% for-loop rather than across the columns, rows and layers in an effort to
% hasten the iteration process, but it actually turned out to take
% marginally longer overall. If you can think of methods to speed up the
% process, please feel free to implement them and make mention of them to
% me.

%% Plotting the Computed Burgers Vectors
% Note that everything has been scaled up by 1e09, e.g. 500e-9 becomes 500,
% to match with how the 3D material shape was coded.

% Plot quivers depicting the coordinate axes. Change the position and size
% to suit your own material specimen.
quiver3(-490,-490,-490,500./0.9,0,0, 'r', 'LineWidth', 1.5);    % x-axis
quiver3(-490,-490,-490,0,500./0.9,0, 'g', 'LineWidth', 1.5);    % y-axis
quiver3(-490,-490,-490,0,0,500./0.9, 'b', 'LineWidth', 1.5);    % z-axis


% Produce meshgrids across which to plot the set of computed Burgers vectors
[mmesh,nmesh,omesh] = meshgrid(m,n,o);

% Plot the computed Burgers vectors. These have been scaled up by an
% additional 10 for improved visibility.
quiver3(mmesh*10^9,nmesh*10^9,omesh*10^9,b_computed(1:numel(n),1:numel(m),1:numel(o),1)*10^10,b_computed(1:numel(n),1:numel(m),1:numel(o),2)*10^10,b_computed(1:numel(n),1:numel(m),1:numel(o),3)*10^10, 'k');


%% Plotting the 3D Material Shape and Dislocation Lines
% The code and data in this section was provided by Felix Hofmann, 
% reference to the associated paper is as follows:
% F. Hofmann et al. “Nanoscale imaging of the full strain tensor of specific dislocations extracted from a bulk sample”. In: Phys. Rev. Materials4.1 (2020), p. 013801. doi:https://doi.org/10.1103/PhysRevMaterials.4.013801.


% To do the same with your own material specimen, you will have to undergo
% a similar method but with your own data.


% load node data:

%plot the dislocation positions...
load('segmentation190123.mat'); % variable in this is called tmp1
tmp1 = segmentation;
nn = 1; % reduced resolution along dislocation line...

% extract dislocation node coordinates and add points outside object where needed. 
ddd = tmp1(1,4); d1xyzt = cell2mat([ddd{:}]);  d1_vect = (d1xyzt(end,:) - d1xyzt(1,:)).*10; 
                        d1xyz = [d1xyzt(1,:)-d1_vect ; d1xyzt(1:nn:end-1,:); d1xyzt(end,:) ;d1xyzt(end,:)+d1_vect];
                        
ddd = tmp1(2,4); d2xyzt = cell2mat([ddd{:}]);  d2_vect = (d2xyzt(end,:) - d2xyzt(1,:)).*10; 
                        d2xyz = [d2xyzt(1,:)-d2_vect ; d2xyzt(1:nn:end-1,:); d2xyzt(end,:) ; d2xyzt(end,:)+d2_vect];
                        
ddd = tmp1(3,4); d3xyzt = cell2mat([ddd{:}]);  d3xyz = [d3xyzt(1,:) + [0 -1000 0]; d3xyzt(1:nn:end-1,:) ; d3xyzt(end,:)];

ddd = tmp1(4,4); d4xyzt = cell2mat([ddd{:}]);  d4xyz = [d4xyzt(1:nn:end-1,:); d4xyzt(end,:); [0 -1000 0] + d4xyzt(end,:)];

ddd = tmp1(5,4); d5xyzt = cell2mat([ddd{:}]);  d5xyz = [d5xyzt(1:nn:end-1,:); d5xyzt(end,:); [0 -1000 0] + d5xyzt(end,:)];

% dislo 3, 4 and 5 don't quite meet in the same place - fix this...
d3xyz(end,:) = d4xyz(1,:); 


% transform dislocation positions from pixel to sample coordiante positions: 
mid_pix = [240./2 240./2 240./2]; % middle pixel position in x, y, z directions 
pix_size = 5*10^-9; %voxel size in m. 
d1xyz_pos = (d1xyz - ones(size(d1xyz,1),1)*mid_pix).*pix_size;
d2xyz_pos = (d2xyz - ones(size(d2xyz,1),1)*mid_pix).*pix_size;
d3xyz_pos = (d3xyz - ones(size(d3xyz,1),1)*mid_pix).*pix_size;
d4xyz_pos = (d4xyz - ones(size(d4xyz,1),1)*mid_pix).*pix_size;
d5xyz_pos = (d5xyz - ones(size(d5xyz,1),1)*mid_pix).*pix_size;

% make sample morphology from reflection amplitude data
for iii = 1: size(Sam_red,1)
amp4d(:,:,:,iii) = Sam_red(iii,1).a_rr;
end
amp_max = squeeze(max(amp4d,[],4));
amp_all = squeeze(mean(amp4d,4));
amp_min = squeeze(min(amp4d,[],4));
clear amp4d
mask_all = amp_all; mask_all(mask_all<0.3) = 0; mask_all(mask_all>0) = 1;


% make rn list of dislocation nodes: 
rn = [d1xyz_pos; d2xyz_pos; d3xyz_pos; d4xyz_pos; d5xyz_pos]+0.1*pix_size; %offset every position by 0.1 od a pixel for numerical stability
rn(:,4) = 7; %make all the nodes fixed. 

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

bv1_dir = -0.5 * [1 1 1]';   bv1 = (bv1_dir(1,1)* n100 + bv1_dir(2,1)* n010 + bv1_dir(3,1)* n001)* a0; bv1n = -bv1./norm(bv1);
bv2_dir = 0.5 * [-1 1 1]';  bv2 = (bv2_dir(1,1)* n100 + bv2_dir(2,1)* n010 + bv2_dir(3,1)* n001)* a0; bv2n = -bv2./norm(bv2);
bv3_dir = -[1 0 0]';         bv3 = (bv3_dir(1,1)* n100 + bv3_dir(2,1)* n010 + bv3_dir(3,1)* n001)* a0; bv3n = -bv3./norm(bv3);
bv4_dir = -0.5 * [1 1 1]';   bv4 = (bv4_dir(1,1)* n100 + bv4_dir(2,1)* n010 + bv4_dir(3,1)* n001)* a0; bv4n = -bv4./norm(bv4);
bv5_dir = -0.5 * [1 -1 -1]'; bv5 = (bv5_dir(1,1)* n100 + bv5_dir(2,1)* n010 + bv5_dir(3,1)* n001)* a0; bv5n = -bv5./norm(bv5);


% work out the segments for the dislocations
Nd1 = size(d1xyz_pos,1);
Nd2 = size(d2xyz_pos,1);
Nd3 = size(d3xyz_pos,1);
Nd4 = size(d4xyz_pos,1);
Nd5 = size(d5xyz_pos,1);

links1(:, 1:2) = [(1:Nd1-1)', (2:Nd1)'];                links1(:,3:5) = ones(Nd1-1,1)*bv1';
links2(:, 1:2) = [(1:Nd2-1)', (2:Nd2)']+Nd1;            links2(:,3:5) = ones(Nd2-1,1)*bv2';
links3(:, 1:2) = [(1:Nd3-1)', (2:Nd3)']+Nd1+Nd2;        links3(:,3:5) = ones(Nd3-1,1)*bv3';
links4(:, 1:2) = [(1:Nd4-1)', (2:Nd4)']+Nd1+Nd2+Nd3;    links4(:,3:5) = ones(Nd4-1,1)*bv4';
links5(:, 1:2) = [(1:Nd5-1)', (2:Nd5)']+Nd1+Nd2+Nd3+Nd4;links5(:,3:5) = ones(Nd5-1,1)*bv5';

links = [links1 ; links2 ; links3; links4; links5];


% Plot dislocation positions overlaid on sample morphology...
amp2plot = amp_all;
pat = patch(isosurface(Xgrid.*10^9, Ygrid.*10^9, Zgrid.*10^9,amp2plot,0.2));
isonormals(Xgrid, Ygrid, Zgrid,amp2plot,pat);
set(pat,'FaceColor','yellow','EdgeColor','none');
alpha(pat,0.1); hold on;
daspect([1,1,1]); view([0 0 1]); axis equal vis3d xy on;
axis([-500 500 -500 500 -500 500]);
camlight;  hold on;

%plot dislocations: 
plot3((d1xyz_pos(2:end-1,1))*1E9, (d1xyz_pos(2:end-1,2))*1E9, (d1xyz_pos(2:end-1,3))*1E9, 'm-', 'LineWidth', 1); hold on
plot3((d2xyz_pos(2:end-1,1))*1E9, (d2xyz_pos(2:end-1,2))*1E9, (d2xyz_pos(2:end-1,3))*1E9, 'c-', 'LineWidth', 1); hold on
plot3((d3xyz_pos(2:end,1))*1E9, (d3xyz_pos(2:end,2))*1E9, (d3xyz_pos(2:end,3))*1E9, 'r-', 'LineWidth', 1); hold on
plot3((d4xyz_pos(1:end-1,1))*1E9, (d4xyz_pos(1:end-1,2))*1E9, (d4xyz_pos(1:end-1,3))*1E9, 'g-', 'LineWidth', 1); hold on
plot3((d5xyz_pos(1:end-1,1))*1E9, (d5xyz_pos(1:end-1,2))*1E9, (d5xyz_pos(1:end-1,3))*1E9, 'b-', 'LineWidth', 1); hold on

%plot burgers vectors (optional):
% bvpl = 120;
% d1pp = (d1xyz_pos(round(end/2),:))*1E9;
% quiver3(d1pp(1,1),d1pp(1,2),d1pp(1,3),bv1n(1,1)*bvpl,bv1n(2,1)*bvpl,bv1n(3,1)*bvpl, 'k', 'LineWidth', 1.5, 'MaxHeadSize',15);
% d2pp = (d2xyz_pos(round(end*0.75),:))*1E9;
% quiver3(d2pp(1,1),d2pp(1,2),d2pp(1,3),bv2n(1,1)*bvpl,bv2n(2,1)*bvpl,bv2n(3,1)*bvpl, 'k', 'LineWidth', 1.5, 'MaxHeadSize',15);
% d3pp = (d3xyz_pos(round(end/2),:))*1E9;
% quiver3(d3pp(1,1),d3pp(1,2),d3pp(1,3),bv3n(1,1)*bvpl,bv3n(2,1)*bvpl,bv3n(3,1)*bvpl, 'k', 'LineWidth', 1.5, 'MaxHeadSize',15);
% d4pp = (d4xyz_pos(round(end/2),:))*1E9;
% quiver3(d4pp(1,1),d4pp(1,2),d4pp(1,3),bv4n(1,1)*bvpl,bv4n(2,1)*bvpl,bv4n(3,1)*bvpl, 'k', 'LineWidth', 1.5, 'MaxHeadSize',15);
% d5pp = (d5xyz_pos(round(end/2),:))*1E9;
% quiver3(d5pp(1,1),d5pp(1,2),d5pp(1,3),bv5n(1,1)*bvpl,bv5n(2,1)*bvpl,bv5n(3,1)*bvpl, 'k', 'LineWidth', 1.5, 'MaxHeadSize',15);
% % add labels to the axes...
% %xlabel('X'); ylabel('Y'); zlabel('Z')
camup([0 1 0]);

axis off


%% Producing a Video Showcasing the Figure
% This section is entirely optional and can be commented out if unwanted.
% This process WILL take a while, and if only the array of computed Burgers
% vectors and the figure are desired then I recommend commenting it out.

% Define the frame-rate of the video:
frame_rate = 30;    % (fps)

% Define the length of the video:
t_length = 10;   % (s)

% Determine the total number of frames. NOTE: There is a strange bug in
% MATLAB where the 'camup' function (used later) does not work beyond half
% or so of the total frames produced. To overcome this, a simple bandaid is
% to double the total number of frames and only use the first half of them.
g = 2*frame_rate*t_length;


% Iterate through each frame:
% (NOTE : the iteration process is timed and the percentage completion is
% also provided)
disp('Generating Frames...')
tic
for f = 1:g
    
    % Define how the POV of the camera changes over time. Change as
    % desired.
    view(-200+0.469*f,-75+0.276*f);
    
    % Set the camera to orientate itself such that the positive y-axis is
    % 'up'. Change as desired.
    camup([0 1 0]);
    
    % This pause is required to ensure each frame is recorded:
    pause(0.0000000000001)
    
    % Store each frame in the movieVector array:
    movieVector(f) = getframe(gcf);
    toc
    
    % Establish percentage completion:
    progress = f/g*100;
    fprintf('%.1f%%\n', progress)
end
disp('Frames Generated')

% Define the name and file-type of the video:
myWriter = VideoWriter('Dislocations','MPEG-4');

myWriter.FrameRate = frame_rate;
% Define the quality of the video:
myWriter.Quality = 95;

% Write the video:
open(myWriter);
writeVideo(myWriter,movieVector(1:int16(g/2)));     % Note only the first half of the frames are used
close(myWriter);