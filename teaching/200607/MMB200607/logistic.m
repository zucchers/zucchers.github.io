% Name:     logistic.m
% Author:   Simone Zuccher
% Created:  16 Apr 2007
% Purpose:  Compute x(k+1) = A*x(k)*(1-x(k))
% Input:    A, number of total iterations and x(1)
% Output:   Plot of x(k) versus k and x(k+1) versus x(k) 
% Modified: 

% Clear all variables
clear

% Change format to long exponential
format long e

% Input the paramter A
A=input('Input A>0 (parameter of x(k+1)=A*x(k)*(1-x(k)) ): ');
if(A<0)
  A=abs(A);
endif

% Input the number of total iterations
m=input('Input N (number of iterations): ');

% Input the initial value. If negative, it will be a random number
s(1)=input('Input 0<=s(1)<=1 (if <0 or >1 then random): ');

% Set s(1) random if out of range
if((s(1)>1) || (s(1)<0))
  s(1)=rand(1);
endif

% Display value of s(1)
disp(s(1));

% Assign x and y needed for 2nd plot
x(1)=s(1);
y(1)=0.0;

% Loop on all points
for n=1:1:m-1
  s(n+1)     = A*s(n)*(1 - s(n));
  disp(s(n+1));
  x(2*n)     = s(n);
  y(2*n)     = s(n+1);
  x(2*n+1)   = s(n+1);
  y(2*n+1)   = s(n+1);
end

% Plots s(n) versus n
gset auto
gset key left Left reverse
plot(s,'+-g;s(n);');

% Wait for keypressed
disp('Please press a key to continue...');
pause();

% Plot f(x), x and path
t=linspace(0,1,500);
plot(t,A*t.*(1-t),'-g;A*x*(1-x);',t,t,'-b;x;',x,y,'+-');

disp('Last 8 points:');
disp(s(m-7:m)');
