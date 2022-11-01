% Simple Upwind scheme

function [Y,x_bins, x] = Simple_UPW(h, epsi, alpha, beta,b, c,  f)
    %h - grid size
    %epsilon - perturbation parameter
    % alpha is left boundary value, beta is right boundary value 
    %b, c, f are the coefficient functions for the problem - usually b, c
    %are constants
    
    x= (0:h:1)';

    x_bins = length(x);

    f_vec = zeros(x_bins, 1);
    f_vec(1) = alpha;    % set left boundary conditions
    f_vec(x_bins) = beta;   % Set right boundary conditons 
    for i=2:1:x_bins-1
        f_vec(i) = f(x(i));
    end

    A = sparse(x_bins, x_bins);
    A(1,1) = 1;
    A(x_bins, x_bins) = 1;

    for j = 2:1:x_bins-1
       A(j,[j-1, j, j+1]) = [(-epsi/(h^2))-(1/(h)*max(b((j-1)*h),0)),c((j-1)*h)+(2*epsi)/(h^2) + (1/h)*abs(b((j-1)*h)),(-epsi/(h^2))+min(b((j-1)*h),0)*(1/(h))];
    end
    
    Y = A\f_vec;
end