## Copyright (C) 2010 Simone Zuccher
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

##
##    function [U, V] = blasius(y,Uinf,nu,x)
##
##    INPUT
##        y:    vector of wall-normal distance on a generally uneven grid
##        Uinf: external streamwise velocity
##        nu:   kinematic viscosity
##        x:    distance from the leading edge,
##    OUTPUT
##        U(y)  stream-wise velocity
##        V(y)  wall-normal velocity
##
##    blasius.m solves Blasius' equation
##        f f'' + 2f''' = 0
##    with BC
##        f(0) = f'(0) = 0 and f'(ymax) = 1
##    by introducing u = f', so as to recast the third-order ODE into system
##          _
##         |
##         |  f u' + 2 u'' = 0
##        <
##         |  f' - u    = 0
##         |_        
##   
##    with BC
##        f(0) = u(0) = 0 and u(ymax) = 1.
##    The last equation for f at ymax is 
##        f'(ymax) - u(ymax) = 0.
##    If you experience convergence problems, make sure to use grid nodes 
##    clustered close to y=0.
##    Uinf=1, nu=1 and x=1 return U and V as functions of eta (similarity
##    variable.

##    Author:   Simone Zuccher <zuccher@sci.univr.it>
##    Created:  18 June 2010
##    Last modified: 


function [U, V] = blasius(y,Uinf,nu,x)
%--------------------------------- Preamble -----------------------------------%
N=length(y); y=y(:); eta=y*sqrt(Uinf/(nu*x));
% Sanity-check on eta
if (min(eta)<0|max(eta)<8|x<=0)
   disp(['error: blasius.m: eta negative or max(eta) < 8 or x<=0. '...
         'Exiting blasius.m...']);
   U=[]; V=[];return
end
%-------------------------------- Constants -----------------------------------%
% These constants are use in the FD code
nv=2; ne=2; uv=0; fv=1; e1=0; e2=1;
% A temporary rectangular matrix which will become a sparse diagonal matrix
smin=-1; smax=1; sizeJac=[N*ne,nv*(smax-smin+2)];
% md: column of rectangular matrix corresponding to the main diagonal
md = -(smin-1)*nv;
%--------------------------------- Newton IC ----------------------------------%
% Initial guess for Newton: use U'(0) = 0.332 and solve a Cauchy problem
F=zeros(size(eta)); U=zeros(size(eta));
j=(2:N-1);
d10(j) = 1./(eta(j+1)-eta(j-1));       d10=d10(:);
d2p(j) = 2*d10(j)./(eta(j+1)-eta(j));  d2p=d2p(:);
d2m(j) = 2*d10(j)./(eta(j)-eta(j-1));  d2m=d2m(:);
d20(j) = -d2p(j)-d2m(j);               d20=d20(:);
U(2)=0.332*eta(2);  % key choice for guess
for j=2:N-1
   F(j)=F(j-1)+(U(j-1)+U(j))/2.*(eta(j)-eta(j-1));
   U(j+1)=-(2*d20(j).*U(j)+(2*d2m(j)-d10(j).*F(j)).*U(j-1))./ ...
          (2*d2p(j)+d10(j).*F(j));
end
F(N)=F(N-1)+(U(N-1)+U(N))*(eta(N)-eta(N-1))/2;
f=zeros(2*N,1); j=(1:N); eq=(j-1)*ne+1; f(eq+uv)=U; f(eq+fv)=F;
j=(1:N); eq=(j-1)*ne+1; f(eq+uv)=U; f(eq+fv)=F;
%-------------------------------- Newton loop ---------------------------------%
err=1;
while(err>1e-10)
   % Set everything back to zero
   Jac=zeros(sizeJac);
   rhs=zeros(sizeJac(1),1);
   %
   % Boundary conditions at y=0
   j=1; eq = (j-1)*ne+1;
   rhs(eq+e1)=f(eq+uv)-0;
   rhs(eq+e2)=f(eq+fv)-0;
   Jac(eq+uv,md-e1+uv)=1;
   Jac(eq+fv,md-e2+fv)=1;
   %
   % Internal points from 2 to N-1
   j=(2:N-1); eq = (j-1)*ne+1;
   rhs(eq+e1)=f(eq+fv).*d10(j).*(f(eq+uv+nv)-f(eq+uv-nv)) + ...
               2*(f(eq+uv+nv).*d2p(j) + f(eq+uv).*d20(j) + f(eq+uv-nv).*d2m(j));
   rhs(eq+e2)=d10(j).*(f(eq+fv+nv)-f(eq+fv-nv)) - f(eq+uv);
   Jac(eq+uv-nv,md-e1+uv-nv)= -f(eq+fv).*d10(j) + 2.*d2m(j);
   Jac(eq+uv,md-e1+uv)= 2*d20(j);
   Jac(eq+uv+nv,md-e1+uv+nv)= f(eq+fv).*d10(j) + 2.*d2p(j);
   Jac(eq+fv,md-e1+fv)= d10(j).*(f(eq+uv+nv)-f(eq+uv-nv));
   Jac(eq+fv-nv,md-e2+fv-nv)= -d10(j);
   Jac(eq+uv,md-e2+uv)= -1;
   Jac(eq+fv+nv,md-e2+fv+nv)= d10(j);
   %
   % Boundary conditions at y=ymax
   j=N; eq = (j-1)*ne+1;
   rhs(eq+e1)=f(eq+uv)-1;
   rhs(eq+e2)=(f(eq+fv)-f(eq+fv-nv))/(eta(N) - eta(N-1))-f(eq+uv);
   Jac(eq+uv,md-e1+uv)=1;
   Jac(eq+fv,md-e2+fv)=1/(eta(N) - eta(N-1));
   Jac(eq+uv,md-e2+uv)=-1;
   Jac(eq+fv-nv,md-e2+fv-nv)=-1/(eta(N) - eta(N-1));
   %
   % Make the Jacobian a sparse diagonal matrix from a rectangular one and solve
   Jac=spdiags(Jac,(1:sizeJac(2))-md,sizeJac(1),sizeJac(1));
   ftmp = -Jac\rhs;
   % Update solution
   f = ftmp + f;
   % Check the absolute error
   err=norm(ftmp);
endwhile
%---------------------------------- Output ------------------------------------%
j=(1:N); eq=(j-1)*ne+1; 
U=Uinf*f(eq+uv); V=Uinf/2*((y/x).*U-sqrt(nu/(Uinf*x))*f(eq+fv));
endfunction
