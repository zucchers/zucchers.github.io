% Name:     esempio4.m
% Author:   Simone Zuccher
% Created:  04 Apr 2007
% Purpose:  Compute the solution of x^k = cos (x/k)
% Input:    Number of total iterations and x0 for nonlinear solver
% Output:   Plot  x^k and cos (x/k)
% Modified: 

% Clear all variables
clear

% Change format to long exponential
format long e

% Input the number of total iterations
m=input('Input N (number of iterations): ');

% Input the initial value for nonlinear solver. 
x0=input('Input x0 (initial value for nonlinear solver): ');

% Create a vector needed for plots
t=linspace(-1.1,1.1,500);

% Loop on all functions
for n=1:1:m
  % set range for plot
  gset xrange[-1.1:1.1]

  % set attributes for plots
  attr1=['-g;t^' int2str(n) ';'];
  attr2=['-b;cos(t/' int2str(n) ');'];
  
  % plot of x^n and cos(x/n)
  plot(t,t.^n,attr1,t,cos(t/n),attr2);
  
  % define the function 
  fun= ["x^" int2str(n) "-cos(x/" int2str(n) ")"];
  
  % compute zero closest to x0 and diplay it
  s(n)=fsolve(fun,x0);
  disp(s(n))
  
  % Wait for keypressed
  disp('Please press a key to continue...');
  pause();
end
