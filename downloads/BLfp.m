## Copyright (C) 2009 Simone Zuccher
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

## Solve the Boundary-layer equations past a flat plate
##
##  |
##  |  u_x + v_y = 0
## <
##  |  u u_x + v u_y = nu u_{yy}
##  |
##
##  with BC:  
##  @ y = 0:     u(x,0) = v(x,0) = 0
##  @ y = inf:   v_y(x,inf) = 0
##  and IC:
##  @ x = 0:     u(0,y) = 1
##  @ x = 0:     v(0,y) = 0
##
##  In the conservative form
##
##  |
##  |  u_x + v_y = 0
## <
##  |  (u^2)_x + (v u)_y = nu u_{yy}
##  |

##
## Since the system is parabolic-in-space, a marching technique is employed
## from x=0 to x=1.


## Author:   Simone Zuccher <zuccher@sci.univr.it>
## Created:  22 May 2009
## Modified: 


clear all
close all
N = 100;     % number of y points
ymax = 10;   % maximum value of y (inf)
xmin = 0;    % initial x
xmax = 1;    % final x-location
nx = 100;     % number of grid points in x
nu = 1;      % viscosity

% Generate y grid clustering points close to the wall
y=((0:(N-1))./(N-1)).^3*ymax;
% Generate x grid clustering points close to x=0
x=xmin + ([0:(nx-1)]./(nx-1)).^2*(xmax-xmin);

% Initial guess for Newton:
% u is 1 except at the wall
% v is zero everywhere
f=zeros(2*N,1);
f(3:2:2*N) = ones(N-1,1);
% Now perturb a little bit to make sure nothing is zero
f = f + 0.01*(0.5-rand(size(f)));

% Solution at previous x
fold = f;

% These constants are use in the FD code
nvar=2;
neq=2;
uv=0;
vv=1;
eq1=0;
eq2=1;


% Loop on x
for ix=2:nx
dx = x(ix) - x(ix-1);

% Set the error to start the Newton's loop
err=1;
while(err>1e-6)
%  Set back everything to zero
   J=spalloc(2*N,2*N,5*2*N);
   rhs=zeros(2*N,1);
   eq=1;
   rhs(eq+eq1)=f(eq+uv)-0;
   rhs(eq+eq2)=f(eq+vv)-0;
   J(eq+eq1,eq+uv)=1;
   J(eq+eq2,eq+vv)=1;
   for i=2:N-1
      d10 = 1/(y(i+1)-y(i-1));
      d1p = 1/(y(i+1)-y(i));
      d1m = 1/(y(i)-y(i-1));
      d2p = 2*d10/(y(i+1)-y(i));
      d20 = -2*d10*(1/(y(i+1)-y(i)) + 1/(y(i)-y(i-1)));
      d2m = 2*d10/(y(i)-y(i-1));
      eq = i*2-1;
      rhs(eq+eq1)=(f(eq+uv)+f(eq+uv-nvar))/2/dx - ...
                  (fold(eq+uv)+fold(eq+uv-nvar))/2/dx + ...
		  d1m*(f(eq+vv)-f(eq+vv-nvar));
      rhs(eq+eq2)=(f(eq+uv)^2 - fold(eq+uv)^2)/dx + ...
                  d10*(f(eq+uv+nvar)*f(eq+vv+nvar) - ...
		       f(eq+uv-nvar)*f(eq+vv-nvar)) ...
		  - nu*(f(eq+uv+nvar)*d2p + f(eq+uv)*d20 + f(eq+uv-nvar)*d2m);
      %
      J(eq+eq1,eq+uv-nvar)= 1/2/dx;
      J(eq+eq1,eq+uv)= 1/2/dx;
      J(eq+eq1,eq+vv-nvar)= -d1m;
      J(eq+eq1,eq+vv)= d1m;
      %
      J(eq+eq2,eq+uv-nvar)= - d10*f(eq+vv-nvar) - nu*d2m;
      J(eq+eq2,eq+uv)= 2*f(eq+uv)/dx - nu*d20;
      J(eq+eq2,eq+uv+nvar)= d10*f(eq+vv+nvar) - nu*d2p;
      J(eq+eq2,eq+vv-nvar)= - d10*f(eq+uv-nvar);
      J(eq+eq2,eq+vv+nvar)= d10*f(eq+uv+nvar);
   end
   eq = N*2-1;
   rhs(eq+eq1)=f(eq+uv)-1;
   rhs(eq+eq2)=(f(eq+vv)-f(eq+vv-nvar))/(y(N) - y(N-1))-0;
   J(eq+eq1,eq+uv)=1;
   J(eq+eq2,eq+vv)=1/(y(N) - y(N-1));
   J(eq+eq2,eq+vv-nvar)=-1/(y(N) - y(N-1));

   df = -J\rhs;
%  Update solution
   f = df + f;
%  Check the error
   err=norm(df);
endwhile
err;
fold = f;
end

UBL=f(1:2:end);
VBL=f(2:2:end);
yBL=y;
plot(UBL,yBL,'-',VBL,yBL,'-')
ylabel('y');
xlabel('u, v');
legend('u','v')
