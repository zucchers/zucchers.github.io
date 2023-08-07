% Name:     esempio3.m
% Author:   Simone Zuccher
% Created:  03 Apr 2007
% Purpose:  Compute x(k+1) = integrate(exp(-t^2),t,0,x(k))
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
s(1)=1.0;

% Display value of s(1)
disp(s(1));

% Assign x and y needed for 2nd plot
x(1)=s(1);
y(1)=0.0;

% Loop on all points
for n=1:1:m-1
  s(n+1)     = sqrt(pi)*erf(s(n))/2.0;
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

disp(s(m))

% Plot f(x), x and path
t=linspace(0,1.5,500);
plot(t,sqrt(pi)*erf(t)/2.0,'-g;sqrt(pi)*erf(x)/2;',t,t,'-b;x;',x,y,'+-');
