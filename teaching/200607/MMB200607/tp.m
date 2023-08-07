% Name:     tp.m
% Author:   Simone Zuccher
% Created:  24 May 2007
% Purpose:  solve the threshold phenomena system
%           u'= u(1+ 1/(1+(u-2)^2) - v)
%           v'= v(u - (v+1))
%           given u0 and v0
% Input:    see file
% Output:   1. plot of u(t) versus v(t) together with the vector field 
%           2. plot time histories of u(t), v(t)
% Modified: 
%
% The non-negative equilibrium points are the following
%
% [u = 2.353209961160612, v = 1.353209961160612]
%

% Clear all variables
clear all;

% Window ranges
xmin=0;
xmax=5;
ymin=0;
ymax=3;


% Equilibrium points
eq = [ 2.682327816337427 1.682327816337427];
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
  xdot(1) = u*(1+ 1/(1+(u-2)^2) - v);
  xdot(2) = v*(u - (v+1));
endfunction

__gnuplot_set__ nokey

setax=[xmin xmax ymin ymax];
axis(setax)

[X, Y] = meshgrid(xmin:dx:xmax, ymin:dy:ymax);
DX = X.*(1+ 1./(1+(X-2).^2) - Y);
DY = Y.*(X - (Y+1));
L = sqrt(DX.^2 + DY.^2);
mytitle=["Phase portrait. Initial conditions: x0=" num2str(x0(1)) \
         ", y0=" num2str(x0(2))];
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
mytitle=["Time histories. Initial conditions: x0=" num2str(x0(1)) \
         ", y0=" num2str(x0(2))];
__gnuplot_set__ auto
__gnuplot_set__ xlabel 't'
__gnuplot_set__ ylabel 'x(t), y(t)'
title(mytitle)
__gnuplot_set__ key

% Plot time histories
plot( t, x(1,:), '-r;x(t);', t, x(2,:), '-g;y(t);')
