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

## Solve Blasius' equation
##
## f f'' + 2f''' = 0
##
## with BC
##
## f(0) = f'(0) = 0 and f'(ymax) = 1.
##
## By introducing u = f', the equation can be recast in the system
##
##  |
##  |  f u' + 2 u'' = 0
## <
##  |  f' - u    = 0
##  |
## with BC
##
## f(0) = u(0) = 0 and u(ymax) = 1.
##
## The last equation for f at ymax is f'(ymax) - u(ymax) = 0

## Author:   Simone Zuccher <zuccher@sci.univr.it>
## Created:  21 May 2009
## Modified: 


%clear all
%close all
N = 200;     % number of y points
ymax = 10;   % maximum value of y (inf)
Uinf = 1;    % U upstream and for y -> inf
nu = 1;      % for comparison with BLfp
xmax=1;      % for comparison with BLfp

% Generate y grid clustering points close to the wall. 
% This is, in fact, eta
y=([0:(N-1)]./(N-1)).^3*ymax;

% Initial guess for Newton:
% u goes from 0 to 1 linearly
% f is zero everywhere
f=zeros(2*N,1);
f(1:2:2*N) = linspace(0,1,N)';
% Now perturb a little bit to make sure nothing is zero
f = f + 0.01*rand(size(f));

% These constants are use in the FD code
nvar=2;
neq=2;
uv=0;
fv=1;
eq1=0;
eq2=1;

% Set the error to start the Newton's loop
err=1;
while(err>1e-10)
%  Set back everything to zero
   J=spalloc(2*N,2*N,5*2*N);
   rhs=zeros(2*N,1);
   eq=1;
   rhs(eq+eq1)=f(eq+uv)-0;
   rhs(eq+eq2)=f(eq+fv)-0;
   J(eq+eq1,eq+uv)=1;
   J(eq+eq2,eq+fv)=1;
   for i=2:N-1
      d10 = 1/(y(i+1)-y(i-1));
      d2p = 2*d10/(y(i+1)-y(i));
      d20 = -2*d10*(1/(y(i+1)-y(i)) + 1/(y(i)-y(i-1)));
      d2m = 2*d10/(y(i)-y(i-1));
      eq = i*2-1;
      rhs(eq+eq1)=f(eq+fv)*d10*(f(eq+uv+nvar)-f(eq+uv-nvar)) + ...
             2*(f(eq+uv+nvar)*d2p + f(eq+uv)*d20 + f(eq+uv-nvar)*d2m);
      rhs(eq+eq2)=d10*(f(eq+fv+nvar)-f(eq+fv-nvar)) - f(eq+uv);
      %
      J(eq+eq1,eq+uv-nvar)= -f(eq+fv)*d10 + 2*d2m;
      J(eq+eq1,eq+uv)= 2*d20;
      J(eq+eq1,eq+uv+nvar)= f(eq+fv)*d10 + 2*d2p;
      J(eq+eq1,eq+fv)= d10*(f(eq+uv+nvar)-f(eq+uv-nvar));
      %
      J(eq+eq2,eq+fv-nvar)= -d10;
      J(eq+eq2,eq+uv)= -1;
      J(eq+eq2,eq+fv+nvar)= d10;
   end
   eq = N*2-1;
   rhs(eq+eq1)=f(eq+uv)-1;
   rhs(eq+eq2)=(f(eq+fv)-f(eq+fv-nvar))/(y(N) - y(N-1))-f(eq+uv);
   J(eq+eq1,eq+uv)=1;
   J(eq+eq2,eq+fv)=1/(y(N) - y(N-1));
   J(eq+eq2,eq+uv)=-1;
   J(eq+eq2,eq+fv-nvar)=-1/(y(N) - y(N-1));

   ftmp = -J\rhs;
%  Upade solution
   f = ftmp + f;
%  Check the error
   err=norm(ftmp);
endwhile

% Change variables according to similarity solution
y = y/sqrt(Uinf/nu/xmax);
U = Uinf*f(1:2:end);
V = 1/sqrt(Uinf*nu/xmax)/2*(y'.*f(1:2:end)-f(2:2:end));
% Plot results
plot(U,y,'-',V,y,'-')
ylabel('\eta');
xlabel('u, v');
legend('u','v')
axis([0 1 0 10])
