## Copyright (C) 2008 Simone Zuccher
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##

## x = odeplot (f, t0, x0, tmin, tmax, dt)
##   Solve the differential equation
## 
##           dx
##           -- = f(x, t)
##           dt
## 
##   with
## 
##           x(t_0) = x_0
##
##   on the interval J=[tmin;tmax]. Vectors t and x are returned as outputs.
##
## Inputs
##   f     Function f(x)
##   t0    Time at which the initial condition x0 is provided 
##   x0    Initial condition at time t0
##   tmin  Lower bound of inteval J
##   tmax  Upper bound of inteval J
##   dt    Time step for numerical integration (delta t)
##
## Example
##
##   [t,x]=odeplot("1/(sqrt(x^2+t^2)-1)",0,0,-.41,.41,.005);
##

## Author:   Simone Zuccher <zuccher@sci.univr.it>
## Created:  08 May 2008
## Modified: 


function [t, x] = odeplot(f, t0, x0, tmin, tmax, dt);

% Sanity checks
if nargin < 6
   usage ("x = odeplot (f, t0, x0, tmin, tmax, dt)");
   return
endif
if tmin > t0
   tmin=t0
endif
if tmax < t0
   tmax=t0
endif

global func
func=f;

% Close all possible plots and change format to long exponential
close all

% Forward and backward time vectors
tf=(t0:dt:tmax);
tb=(t0:-dt:tmin);

% Forward and backward solutions
xf = lsode("f", x0, tf)';
xb = lsode("f", x0, tb)';
key=[";f(x) = " f ", t0 = " num2str(t0) ", x0 = " num2str(x0) ";"];

% Plot
plot(tf,xf,key,tb,xb,"-r;;");

% Return complete vectors of time and solution
x=[xb((end:-1:1)) xf];
t=[tb((end:-1:1)) tf];
return
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the ODE function f(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xdot=f(x, t)
  global func;
  eval(["xdot=" func ";"]);
endfunction
