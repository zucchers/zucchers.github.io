% Name:     esempio2.m
% Author:   Simone Zuccher
% Created:  03 Apr 2007
% Purpose:  Compute x(k+1) = max(1/4,x(k)^2)
% Input:    Number of total iterations and x(1)
% Output:   Plot of x(k) versus k and x(k+1) versus x(k) 
% Modified: 

% Clear all variables
clear

% Change format to long exponential
format long e

% Input the number of total iterations
m=input('Input N (number of iterations): ');

% Input the initial value. 
s(1)=input('Input s(1): ');

% Display value of s(1)
disp(s(1));

% Assign x and y needed for 2nd plot
x(1)=s(1);
y(1)=0.0;

% Loop on all points
for n=1:1:m-1
  s(n+1)     = max(1/4,s(n)^2);
  disp(s(n+1));
  x(2*n)     = s(n);
  y(2*n)     = s(n+1);
  x(2*n+1)   = s(n+1);
  y(2*n+1)   = s(n+1);
end

% Plots s(n) versus n
gset auto
plot(s,'+-g;s(n);');

% Wait for keypressed
disp('Please press a key to continue...');
pause();

% Display value of s(m)
disp(s(m));

% Plot f(x), x and path
t=linspace(-1.5,1.5,500);
plot(t,max(1/4,t.^2),'-g;max(1/4,t^2);',t,t,'-b;x;',x,y,'+-');
