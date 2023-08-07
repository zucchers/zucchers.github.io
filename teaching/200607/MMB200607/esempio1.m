% Name:     esempio1.m
% Author:   Simone Zuccher
% Created:  03 Apr 2007
% Purpose:  Compute x(k+1) = x(k) + sin(x(k))
% Input:    Number of total iterations and x(1)
% Output:   Plot of x(k) versus k and x(k+1) versus x(k) 
% Modified: 

% Clear all variables
clear

% Change format to long exponential
format long e

% Input the number of total iterations
m=input('Input N (number of iterations): ');

% Input the initial value. If negative, it will be a random number
s(1)=input('Input 0<s(1)<pi (if <0 or >pi then random): ');

% Set s(1) random if out of range
if((s(1)>pi) || (s(1)<0))
  s(1)=pi*rand(1);
endif

% Display value of s(1)
disp(s(1));

% Assign x and y needed for 2nd plot
x(1)=s(1);
y(1)=0.0;

% Loop on all points
for n=1:1:m-1
  s(n+1)     = s(n)+sin(s(n));
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

% Plot f(x), x and path
t=linspace(0,pi,500);
plot(t,t+sin(t),'-g;x+sin(x);',t,t,'-b;x;',x,y,'+-');
