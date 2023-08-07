% Name:     esempio7.m
% Author:   Simone Zuccher
% Created:  03 Apr 2007
% Purpose:  Compute x(k+1) = x(k) + [x(k-1)]^2
% Input:    Number of total iterations and x(2)
% Output:   Plot of x(k) versus k and x(k+1) versus x(k) 
% Modified: 

% Clear all variables
clear

% Change format to long exponential
format long e

% Input the number of total iterations
m=input('Input N (number of iterations): ');

% Input the initial value. 
s(1)=0.;
s(2)=input('Input s(2) (second initial value): ');

% Display value of s(1) and s(2)
disp(s(1));
disp(s(2));

% Loop on all points
for n=2:1:m-1
  s(n+1)     = s(n) + s(n-1)^2;
  disp(s(n+1));
end

% Plots s(n) versus n
gset auto
plot(s,'+-g;s(n);');

disp('');
disp(s(m))
