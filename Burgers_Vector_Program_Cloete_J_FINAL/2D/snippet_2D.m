function strains_snippet = snippet_2D(strains,xlims,ylims,xmin,ymax,interval)
% snippet_2D
% Extracts a snippet of the total 4D displacement gradient array for use as
% the input to the find_burgers_vector function

% Note that the most negative x-coordinates and the most positive
% y-coordinate for the sample data are also required

% Find the corresponding element positions for the x- and y-limits of the
% snippet:
x_start = int16((xlims(1) - xmin)/interval + 1);
% Note that we use ymax and subtract ylims, because the first element in
% the y-direction has been specified to correspond to the most positive
% y-coordinate.
y_start = int16((ymax - ylims(1))/interval + 1);

x_end = int16((xlims(2) - xmin)/interval + 1);
y_end = int16((ymax - ylims(2))/interval + 1);


% Extract the snippet from the total 4D displacement gradient array:
strains_snippet = strains(y_start:y_end,x_start:x_end,:,:);