% Name:     population.m
% Author:   Simone Zuccher
% Created:  17 Apr 2007
% Purpose:  Compute the logistic map
%           x(k+1) = A*x(k)*(1-x(k))
%           computed up to n  given A and  x_0
% Input:    A, , x_0 and number of total iterations n
% Output:   Array of x values 
% Modified: 

function x =  population(A, x0, n)
  x = zeros(n, 1);
  x(1) = x0;
  for k = 1:n-1
    x(k + 1) = A * x(k) * (1 - x(k));
  end
