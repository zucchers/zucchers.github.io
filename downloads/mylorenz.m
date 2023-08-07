function [t, x] = mylorenz(dt,tmax,param,x0,rep);
% This version has been tested on:
% + GNU Octave 3.0
%
% mylorenz.m solves the classic Lorenz system
%      _
%     |
%     |  d x
%     |  --- =  a ( y - x )
%     |  d t
%     |
%     |  d y
%    <   --- =  - x z + b x - y
%     |  d t
%     |
%     |  d z
%     |  --- =  x y - c z 
%     |  d t
%     |_        
%
% given the initial condition [x0 y0 w0] and the parameters [a b c], and plots:
%      (a) the trajectory (x,y,z) in the phase portrait 
%      (b) the time histories of x(t), y(t), z(t).
%
%
% Inputs
%   dt    Time-step for time integration (can be negative for backward-in-time
%         integration)
%   tmax  Maximum time (must be negative if dt < 0 )
%   param Vector of parameters [a b c]
%   x0    Vector of initial conditions x(0)= [x0 y0 z0]
%   rep   Flag: 1 - replot the new trajectory on the existing plot
%               0 - plot everything from scratch
%
%
% Outputs
%   t     Vector of time
%   x     Solution, (x in 1st column, y in the 2nd, z in the 3rd)
%
% Examples
%
% 1. Plot the vector field, equilibria and a trajectory from scratch:
%
%    [t,x] = mylorenz(0.01,10,[10 3/2 8/3],[-1 .5 1],0);
%
% 2. Plot a new trajectory on the existing plot:
%
%    [t,x] = mylorenz(0.01,10,[10 3/2 8/3],[0 0 -1],1);
%
%

% Copyright (C) 2008 Simone Zuccher
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% Author:   Simone Zuccher <zuccher@sci.univr.it>
% Created:  10 Jun 2008
% Modified: 

% Sanity check
if nargin < 5
   usage ('[t, X] = mylorenz(x0,dt,tmax,aa,rep);');
   return
end

if(rep == 0) 
   close all; 
end


%disp('Initial condition:');
%disp(x0);

% Model constant
global aa;
aa=param;

bstar2=aa(1)*(aa(1)+aa(3)+3)/(aa(1)-aa(3)-1);
%disp('b**:')
%disp(bstar2)

% Equilibrium points
eq = [ 0 0 0];
eq = [eq; sqrt(aa(2)*aa(3)-aa(3)) sqrt(aa(2)-1)*sqrt(aa(3)) aa(2)-1];
eq = [eq; -sqrt(aa(2)*aa(3)-aa(3)) -sqrt(aa(2)-1)*sqrt(aa(3)) aa(2)-1];

%disp('Equilibrium points:');
%disp(eq);

% Time parameters
tmin=0;
% Create time 
t = tmin:dt:tmax;
t=t';


if (exist('lsode')>0)
% Octave
   t = 0:dt:tmax;
   x =  lsode('dsys', x0, t);
elseif(exist('ode45')>0)
% Matlab
   [t,x]=ode45('dsys',[0 tmax],x0');
end




figure(1)
if(rep == 1)
   hold on
end   
plot3(x(:,1),x(:,2),x(:,3),'-r',x0(1),x0(2),x0(3),'or');
if(~(abs(imag(eq))>0))
   hold on
   plot3(eq(:,1),eq(:,2),eq(:,3),'ob');
else
   disp('WARNING: origin is the only physical equilibrium point')
end
     
view(60,30)
xlabel('x');ylabel('y');zlabel('z');
hold off

figure(2)
plot(t,x(:,1),t,x(:,2),t,x(:,3))
xlabel('t');ylabel('x(t), y(t), z(t)');
legend('x(t)', 'y(t)', 'z(t)')


return




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
