function b_computed = find_burgers_vector_repeated_2D(interval,xlims,ylims,beta,i)
% find_burgers_vector_repeated
% finds a set of Burgers vectors from the 4D displacement gradient array

% The outermost Burgers circuit is as defined in the program's associated
% paper, with the exception of being reduced to a 2D circuit in the x-y
% plane.

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

% Create an empty column vector in which to store the Burgers vectors:
b_computed = zeros(3,(i+1));
for m = 0:1:i
    % Perform the numerical integration around the Burgers circuit:
    % x-component:
    b_computed(1,m+1) = -simps(xvector((1+m):(end-m)),beta((1+m),(1+m):(end-m),1,1))...
        + -simps(yvector((1+m):(end-m)),beta((1+m):(end-m),(1+m),1,2))...
        + simps(xvector((1+m):(end-m)),beta((end-m),(1+m):(end-m),1,1))...
        + simps(yvector((1+m):(end-m)),beta((1+m):(end-m),(end-m),1,2));
    
    % y-component:
    b_computed(2,m+1) = -simps(xvector((1+m):(end-m)),beta((1+m),(1+m):(end-m),2,1))...
        + -simps(yvector((1+m):(end-m)),beta((1+m):(end-m),(1+m),2,2))...
        + simps(xvector((1+m):(end-m)),beta((end-m),(1+m):(end-m),2,1))...
        + simps(yvector((1+m):(end-m)),beta((1+m):(end-m),(end-m),2,2));
    
    % z-component:
    b_computed(3,m+1) = -simps(xvector((1+m):(end-m)),beta((1+m),(1+m):(end-m),3,1))...
        + -simps(yvector((1+m):(end-m)),beta((1+m):(end-m),(1+m),3,2))...
        + simps(xvector((1+m):(end-m)),beta((end-m),(1+m):(end-m),3,1))...
        + simps(yvector((1+m):(end-m)),beta((1+m):(end-m),(end-m),3,2));
end