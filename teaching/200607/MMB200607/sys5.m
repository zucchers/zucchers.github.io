% Name:     sys5.m
% Author:   Simone Zuccher
% Created:  03 May 2007
% Purpose:  Solve the ODE system
%           x'=y
%           y'=(1 - x^2 - y^2)y - x
%           given the initial condition x0 and y0
% Input:    Initial condition [x0 y0], final time tfinal
% Output:   1. plot of x(t) versus y(t) together with the vector field 
%           2. plot time histories of x(t) and y(t)
% Modified: 
%
% The equilibrium points are the following:
%
% [x = 0, y = 0], 
%

% Clear all variables
clear all;

% Window ranges
xmin=-2;
xmax=2;
ymin=-2;
ymax=2;

% Equilibrium points
eq = [0 0];
disp('Equilibrium points:');
disp(eq);

% Set initial conditions
x0=input('Insert initial conditions [x0 y0]: ');

% Set final time for integration
tmax=input('Insert final time: ');


disp('Initial condition:');
disp(x0);

% Time parameters
tmin=0;
dt=.01;
% Create time 
t = tmin:dt:tmax;

% dx and dy used only for vectors
dx=abs(xmax-xmin)/30;
dy=abs(ymax-ymin)/30;
% rescales vector size
scale=0.027*max(abs(xmax-xmin),abs(ymax-ymin));

% Definition of the dynamical system
function xdot=dsys(x, t)
  u = x(1);
  v = x(2);
  xdot(1) = v;
  xdot(2) = (1-u^2-v^2)*v - u;
endfunction

__gnuplot_set__ nokey
axis([xmin xmax ymin ymax]);

[X, Y] = meshgrid(xmin:dx:xmax, ymin:dy:ymax);

DX = Y;
DY = (1-X.^2 - Y.^2).*Y - X;
L = sqrt(DX.^2 + DY.^2);
mytitle=["Phase portrait. Initial conditions: x0=" num2str(x0(1))\
         ", y0=" num2str(x0(2)) "."];
__gnuplot_set__ nokey
__gnuplot_set__ xlabel 'x(t)'
__gnuplot_set__ ylabel 'y(t)'
title(mytitle)

% Plot vector field
quiver(X, Y, scale*DX./L, scale*DY./L)
hold on;

% Plot all equilibrium points
plot( eq(:,1), eq(:,2), '*k')

x =  lsode("dsys", x0, t)';
plot( x(1,1), x(2,1), '*k', x(1,:), x(2,:), '-r')
hold off;

% Wait for keypressed
disp('Please press a key to continue...');
pause();
mytitle=["Time histories. Initial conditions: x0=" num2str(x0(1))\
         ", y0=" num2str(x0(2)) "."];
__gnuplot_set__ auto
__gnuplot_set__ xlabel 't'
__gnuplot_set__ ylabel 'x(t), y(t)'
title(mytitle)
__gnuplot_set__ key

% Plot time histories
plot( t, x(1,:), '-r;x(t);', t, x(2,:), '-g;y(t);')
