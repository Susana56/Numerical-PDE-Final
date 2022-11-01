% Upwind scheme with artificial diffusion

function [Y,x_bins, x] = UPW_AD(h, epsi, alpha, beta,b, c,  f, sigma)
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
       q = (b(x(j))*h)/(2*epsi);
       s = sigma(q);
       A(j,[j-1, j, j+1]) = [(-epsi*s/(h^2))-(b(x(j))/(2*h)),c(x(j))+(2*epsi*s)/(h^2),(-epsi*s/(h^2))+(b(x(j)))/(2*h)];
    end
    
    Y = A\f_vec;
end