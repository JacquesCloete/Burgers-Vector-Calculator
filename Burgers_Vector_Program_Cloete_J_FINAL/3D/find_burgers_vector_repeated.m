function b_computed = find_burgers_vector_repeated(interval,xlims,ylims,zlims,beta,i,j)
% find_burgers_vector_repeated
% finds a set of Burgers vectors from the 5D displacement gradient array

% The outermost Burgers circuit is as defined in the program's associated
% paper.

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
zvector = zlims(1):interval:zlims(2);

% Create an empty column vector in which to store the Burgers vectors:
b_computed = zeros(3,(i+1)*(j+1));
p = 1;
for n = 0:1:j
    for m = 0:1:i
        % Perform the numerical integration around the Burgers circuit:
        % x-component:
        b_computed(1,p) = -simps(xvector((1+n):(end-m)),beta((1+m),(1+n):(end-m),(1+n),1,1))...
            + simps(zvector((1+n):(end-m)),beta((1+m),(1+n),(1+n):(end-m),1,3))...
            + -simps(yvector((1+m):(end-n)),beta((1+m):(end-n),(1+n),(end-m),1,2))...
            + simps(xvector((1+n):(end-m)),beta((end-n),(1+n):(end-m),(end-m),1,1))...
            + -simps(zvector((1+n):(end-m)),beta((end-n),(end-m),(1+n):(end-m),1,3))...
            + simps(yvector((1+m):(end-n)),beta((1+m):(end-n),(end-m),(1+n),1,2));
        
        % y-component:
        b_computed(2,p) = -simps(xvector((1+n):(end-m)),beta((1+m),(1+n):(end-m),(1+n),2,1))...
            + simps(zvector((1+n):(end-m)),beta((1+m),(1+n),(1+n):(end-m),2,3))...
            + -simps(yvector((1+m):(end-n)),beta((1+m):(end-n),(1+n),(end-m),2,2))...
            + simps(xvector((1+n):(end-m)),beta((end-n),(1+n):(end-m),(end-m),2,1))...
            + -simps(zvector((1+n):(end-m)),beta((end-n),(end-m),(1+n):(end-m),2,3))...
            + simps(yvector((1+m):(end-n)),beta((1+m):(end-n),(end-m),(1+n),2,2));
        
        % z-component:
        b_computed(3,p) = -simps(xvector((1+n):(end-m)),beta((1+m),(1+n):(end-m),(1+n),3,1))...
            + simps(zvector((1+n):(end-m)),beta((1+m),(1+n),(1+n):(end-m),3,3))...
            + -simps(yvector((1+m):(end-n)),beta((1+m):(end-n),(1+n),(end-m),3,2))...
            + simps(xvector((1+n):(end-m)),beta((end-n),(1+n):(end-m),(end-m),3,1))...
            + -simps(zvector((1+n):(end-m)),beta((end-n),(end-m),(1+n):(end-m),3,3))...
            + simps(yvector((1+m):(end-n)),beta((1+m):(end-n),(end-m),(1+n),3,2));
        
        p = p+1;
    end
end