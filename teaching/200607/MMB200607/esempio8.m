% Name:     esempio8.m
% Author:   Simone Zuccher
% Created:  12 Apr 2007
% Purpose:  Compute x(k+1) = log( 1 + x(k) )
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
s(1)=input('Input s(1) (initial value): ');

% Set s(1) positive if negative
if(s(1)<0)
  s(1)=abs(s(1));
endif

% Display value of s(1)
disp(s(1));

% Assign x and y needed for 2nd plot
x(1)=s(1);
y(1)=0.0;

% Loop on all points
for n=1:1:m-1
  s(n+1)     = log(1.+s(n));
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
t=linspace(-.5,2,500);
plot(t,log(1.+t),'-g;log(1+t);',t,t,'-b;x;',x,y,'+-');
