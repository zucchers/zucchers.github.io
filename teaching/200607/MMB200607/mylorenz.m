% Name:     mylorenz.m
% Author:   Simone Zuccher
% Created:  24 May 2007
% Purpose:  solve the classic Lorenz system
%           u'= a(v-u)
%           v'= -uw+bu -v
%           z'= uv-cw
%           given u0 and v0
% Input:    see file
% Output:   1. plot of u(t), v(t), w(t) 
%           2. plot time histories of u(t), v(t), w(t)
% Modified: 
%
% The non-negative equilibrium points are the following
%
% [u = 0, v = 0, w = 0]
%
% [u = SQRT(b c - c), v = SQRT(b - 1) SQRT(c), w = b - 1],
%
% [u = - SQRT(b c - c), v = - SQRT(b - 1) SQRT(c), w = b - 1],
%
%

% Clear all variables
clear all;

% Model constant
global aa;

% Set initial conditions
aa=input('Insert coefficients [a b c]: ');
%aa=[10 28 8./3.]

% Equilibrium points
eq = [ 0 0 0];
eq = [eq; sqrt(aa(2)*aa(3)-aa(3)) sqrt(aa(2)-1)*sqrt(aa(3)) aa(2)-1];
eq = [eq; -sqrt(aa(2)*aa(3)-aa(3)) -sqrt(aa(2)-1)*sqrt(aa(3)) aa(2)-1];

disp('Equilibrium points:');
disp(eq);


% Set initial conditions
x0=input('Insert initial conditions [x0 y0]: ');
%x0=[3 15 1];

% Set final time for integration
tmax=input('Insert final time: ');
%tmax=100;

disp('Initial condition:');
disp(x0);

% Time parameters
tmin=0;
dt=.01;
% Create time 
t = tmin:dt:tmax;

% Definition of the dynamical system
function xdot=dsys(x, t)
  global aa;
  u = x(1);
  v = x(2);
  w = x(3);
  xdot(1) = aa(1)*(v-u);
  xdot(2) = -u*w+aa(2)*u - v;
  xdot(3) = u*v-aa(3)*w;
endfunction

x =  lsode("dsys", x0, t');

mytitle=["Phase portrait. Initial conditions: x0=" num2str(x0(1)) \
         ", y0=" num2str(x0(2)) ", z0=" num2str(x0(3)) \
	 ", a=" num2str(aa(1)) ", b=" num2str(aa(2)) \
	 ", c=" num2str(aa(3))];
__gnuplot_set__ nokey
__gnuplot_set__ xlabel 'x(t)'
__gnuplot_set__ ylabel 'y(t)'
__gnuplot_set__ zlabel 'z(t)'
__gnuplot_set__ parametric
__gnuplot_set__ view 
%__gnuplot_set__ view 120, 30
title(mytitle)
gsplot x
