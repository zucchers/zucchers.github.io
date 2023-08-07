% Name:     logisticplot.m
% Author:   Simone Zuccher
% Created:  17 Apr 2007
% Purpose:  Plot last p_max iterations of the logistic map
%           x(k+1) = A*x(k)*(1-x(k))
%           computed up to t_max  for A_min < A < A_max with x_0 = .1
% Input:    A, number of total iterations and x(1)
% Output:   Plot of bifuration diagram
% Modified: 

% Clear all variables
clear

% Change format to long exponential
format long e

% Set A_min and A_max and number of A-values
A_min = 2.8; 
A_max = 4.0;
n = 1000;

% Set maximum number of iterations (we need to stop sooner or later...)
t_max = 1000;

% Set how many iterations from the last are shown for each value of A 
p_max = 100;

% Set the initial point (arbitrary)
x0 = 0.1;

% Create vector of A and matrix for final plot
A = linspace(A_min, A_max, n);
pop = zeros(p_max, n);

% Compute the population for each value of A
for k = 1:n
  x = population(A(k), x0, t_max);
  % Retain only the last p_max iterations
  pop(:, k) = x(t_max-p_max+1:t_max);
end

% Set no key
gset nokey;
% Set the title on x
gset xlabel 'A';
% Set the range of x
gset xrange[2.8:4]
% Generate the plot
plot(A, pop, 'b.');
