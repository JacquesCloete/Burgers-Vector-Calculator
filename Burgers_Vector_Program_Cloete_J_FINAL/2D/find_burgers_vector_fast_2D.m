function b_computed = find_burgers_vector_fast_2D(interval,xlims,ylims,beta)
% find_burgers_vector_fast
% finds the Burgers vector from the 4D displacement gradient array

% The Burgers circuit is as defined in the program's associated paper, with
% the exception of being reduced to a 2D circuit in the x-y plane.

% Note that minus signs have been added to some of the numerical integrals
% below; this is because the trapz function cannot distinguish between the
% processes of integrating from limits a->b and from limits b->a (one is
% the negative of the other)

% To speed up the process when many iterations are required, I have
% implemented a faster version of the trapz function:
% Umberto Picchini (2021). Fast Trapezoidal Integration (https://www.mathworks.com/matlabcentral/fileexchange/8644-fast-trapezoidal-integration), MATLAB Central File Exchange. Retrieved August 17, 2021.

% Define vectors for the variables across which to numerically integrate:
xvector = xlims(1):interval:xlims(2);
yvector = ylims(2):interval:ylims(1);

% Create an empty column vector in which to store the Burgers vector:
b_computed = zeros(3,1);


% Perform the numerical integration around the Burgers circuit:
% x-component:
b_computed(1) = -trapzfu2(xvector,beta(1,:,1,1))...
+ -trapzfu2(yvector,beta(:,1,1,2))...
+ trapzfu2(xvector,beta(end,:,1,1))...
+ trapzfu2(yvector,beta(:,end,1,2));

% y-component:
b_computed(2) = -trapzfu2(xvector,beta(1,:,2,1))...
+ -trapzfu2(yvector,beta(:,1,2,2))...
+ trapzfu2(xvector,beta(end,:,2,1))...
+ trapzfu2(yvector,beta(:,end,2,2));

% z-component:
b_computed(3) = -trapzfu2(xvector,beta(1,:,3,1))...
+ -trapzfu2(yvector,beta(:,1,3,2))...
+ trapzfu2(xvector,beta(end,:,3,1))...
+ trapzfu2(yvector,beta(:,end,3,2));

