% Name:     realLV.m
% Author:   Simone Zuccher
% Created:  24 May 2007
% Purpose:  solve the realistic Lotka-Volterra system
%           u'= u(1-u) - auv/(u+d)
%           v'= bv(1 - v/u)
%           given u0 and v0
% Input:    see file
% Output:   1. plot of u(t) versus v(t) together with the vector field 
%           2. plot time histories of u(t), v(t)
% Modified: 
%
% The non-negative equilibrium points are the following
%
%                            2                  2
%       SQRT(d  + (2 a + 2) d + a  - 2 a + 1) - d - a + 1
% [ u = -------------------------------------------------, 
%                               2
%
%                            2                  2
%       SQRT(d  + (2 a + 2) d + a  - 2 a + 1) - d - a + 1
%   v = ------------------------------------------------- ] 
%                               2
%

% Clear all variables
clear all;

% Window ranges
xmin=0;
xmax=1;
ymin=0;
ymax=.5;

% Model constant
global aa;

% Give instructions
disp('');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
disp('This script solves a real version of the Lotka-Volterra system:');
disp('  u''= u(1-u) - auv/(u+d)');
disp('  v''= bv(1 - v/u)');
disp('given a, b, d ');

% Set initial conditions
aa=input('Insert constants [a b d]: ');

% Compute b_stab (b<b_stab ensures stability)
bstab=(aa(1)-sqrt((1-aa(1)-aa(3))^2+4*aa(3)))*\
      (1+aa(1)+aa(3)-sqrt((1-aa(1)-aa(3))^2+4*aa(3)))/(2*aa(1));
      
if (aa(2)<bstab) 
   disp('WARNING: a, b, d out of the stability region');
endif

% Equilibrium points
eq = [ (sqrt(aa(3)^2+(2*aa(1)+2)*aa(3)+aa(1)^2-2*aa(1)+1)-aa(3)-aa(1)+1)/2\
       (sqrt(aa(3)^2+(2*aa(1)+2)*aa(3)+aa(1)^2-2*aa(1)+1)-aa(3)-aa(1)+1)/2];
%eq = [eq; aa(3)/aa(4) aa(1)/aa(2)];
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
  global aa;
  u = x(1);
  v = x(2);
  xdot(1) = u*(1-u) - aa(1)*u*v/(u+aa(3));
  xdot(2) = aa(2)*v*(1 - v/u);
endfunction

__gnuplot_set__ nokey

setax=[xmin xmax ymin ymax];
axis(setax)

[X, Y] = meshgrid(xmin:dx:xmax, ymin:dy:ymax);
DX = X.*(1-X) - aa(1).*X.*Y./(X+aa(3));
DY = aa(2).*Y.*(1 - Y./X);
L = sqrt(DX.^2 + DY.^2);
mytitle=["Phase portrait. Initial conditions: x0=" num2str(x0(1)) \
         ", y0=" num2str(x0(2)) \
	 "; a=" num2str(aa(1)) ", b=" num2str(aa(2))\
	 ", c=" num2str(aa(3))];
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
         ", y0=" num2str(x0(2)) \
	 "; a=" num2str(aa(1)) ", b=" num2str(aa(2))\
	 ", c=" num2str(aa(3))];
__gnuplot_set__ auto
__gnuplot_set__ xlabel 't'
__gnuplot_set__ ylabel 'x(t), y(t)'
title(mytitle)
__gnuplot_set__ key

% Plot time histories
plot( t, x(1,:), '-r;x(t);', t, x(2,:), '-g;y(t);')
