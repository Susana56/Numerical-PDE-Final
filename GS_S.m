% Gushchin-Shchennikov scheme

function [Y,x_bins, x] = GS_S(h, epsi, alpha, beta,b, c,  f)
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
        f_vec(i) = f(x(i)-h/2);
    end

    A = sparse(x_bins, x_bins);
    A(1,1) = 1;
    A(x_bins, x_bins) = 1;
    A(2, [1,2,3]) = [(-epsi/(h^2))-(b(x(2)-(h/2)))/h + (c(x(2)-(h/2)))/2 , (c(x(2)-(h/2)))/2 + (2*epsi)/(h^2) + (b(x(2)-(h/2)))/h ,(-epsi/(h^2))];
    
    for j = 3:1:x_bins-1
       A(j,[j-2,j-1, j, j+1]) = [(-epsi/(2*h^2)) , (c(x(j)-(h/2)))/2 + (epsi)/(2*h^2) - (b(x(j)-(h/2)))/h ,(epsi/(2*h^2)) + (b(x(j)-(h/2)))/h + (c(x(j)-(h/2)))/2, (-epsi/(2*h^2))];
    end

    Y = A\f_vec;
end 