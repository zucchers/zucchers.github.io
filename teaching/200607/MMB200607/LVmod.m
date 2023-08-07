% Name:     LVmod.m
% Author:   Simone Zuccher
% Created:  17 May 2007
% Purpose:  solve the modified Lotka-Volterra system
%           u'=u(1-u)-kuv^2
%           v'=v(1-v)-kvu^2
%           given u0 and v0, and k>3/4
% Input:    see file
% Output:   1. plot of u(t) versus v(t) together with the vector field 
%           2. plot time histories of u(t), v(t)
% Modified: 
%
% The equilibrium points are the following, but we consider only u0 and v0
% non-negative
%
% [u = 0, v = 0], 
%
% [u = 1, v = 0], 
%
% [u = 0, v = 1],
%
%        sqrt(4 k + 1) + 1        sqrt(4 k + 1) + 1
% [u = - -----------------, v = - -----------------],
%              2 k                      2 k
%
%      sqrt(4 k + 1) - 1      sqrt(4 k + 1) - 1
% [u = -----------------, v = -----------------],
%             2 k                    2 k
%
%      sqrt(4 k - 3) + 1        sqrt(4 k - 3) - 1
% [u = -----------------, v = - -----------------],
%             2 k                      2 k
%
%        sqrt(4 k - 3) - 1      sqrt(4 k - 3) + 1
% [u = - -----------------, v = -----------------]
%               2 k                    2 k
%

% Clear all variables
clear all;

% Window ranges
xmin=0;
xmax=1;
ymin=0;
ymax=1;


% Model constant
global k;
% Set model constant
k=input('Insert k: ');

% Equilibrium points
eq = [0 0];
eq = [eq; 0 1];
eq = [eq; 1 0];
eq= [eq; -(sqrt(4 * k + 1) + 1)/2/k -(sqrt(4 * k + 1) + 1)/2/k];
eq= [eq; (sqrt(4 * k + 1) - 1)/2/k (sqrt(4 * k + 1) - 1)/2/k];
eq= [eq; (sqrt(4 * k - 3) + 1)/2/k -(sqrt(4 * k - 3) - 1)/2/k];
eq= [eq; -(sqrt(4 * k - 3) - 1)/2/k (sqrt(4 * k - 3) + 1)/2/k];
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
  global k;
  u = x(1);
  v = x(2);
  xdot(1) = u*(1-u)-k*u*v^2;
  xdot(2) = v*(1-v)-k*v*u^2;
endfunction

__gnuplot_set__ nokey

setax=[xmin xmax ymin ymax];
axis(setax)

[X, Y] = meshgrid(xmin:dx:xmax, ymin:dy:ymax);

DX = X.*(1-X)-k*X.*(Y.^2) ;
DY = Y.*(1-Y)-k*Y.*(X.^2);
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
