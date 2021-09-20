function b_computed = find_burgers_vector_accurate_2D(interval,xlims,ylims,beta)
% find_burgers_vector_fast
% finds the Burgers vector from the 4D displacement gradient array

% The Burgers circuit is as defined in the program's associated paper, with
% the exception of being reduced to a 2D circuit in the x-y plane.

% Note that minus signs have been added to some of the numerical integrals
% below; this is because the trapz function cannot distinguish between the
% processes of integrating from limits a->b and from limits b->a (one is
% the negative of the other)

% In an attempt to improve the accuracy of the numerical integration, I
% have implemented a version of trapz that uses Simpson's Rule:
% Damien Garcia (2021). Simpson's rule for numerical integration (https://www.mathworks.com/matlabcentral/fileexchange/25754-simpson-s-rule-for-numerical-integration), MATLAB Central File Exchange. Retrieved August 18, 2021.

% Define vectors for the variables across which to numerically integrate:
xvector = xlims(1):interval:xlims(2);
yvector = ylims(2):interval:ylims(1);

% Create an empty column vector in which to store the Burgers vector:
b_computed = zeros(3,1);


% Perform the numerical integration around the Burgers circuit:
% x-component:
b_computed(1) = -simps(xvector,beta(1,:,1,1))...
+ -simps(yvector,beta(:,1,1,2))...
+ simps(xvector,beta(end,:,1,1))...
+ simps(yvector,beta(:,end,1,2));

% y-component:
b_computed(2) = -simps(xvector,beta(1,:,2,1))...
+ -simps(yvector,beta(:,1,2,2))...
+ simps(xvector,beta(end,:,2,1))...
+ simps(yvector,beta(:,end,2,2));

% z-component:
b_computed(3) = -simps(xvector,beta(1,:,3,1))...
+ -simps(yvector,beta(:,1,3,2))...
+ simps(xvector,beta(end,:,3,1))...
+ simps(yvector,beta(:,end,3,2));

