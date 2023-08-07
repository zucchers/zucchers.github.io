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

## x = cobweb (f, N, x1, xmin, xmax)
##   Generate the time history and cobwebbing plot of the one-step map
##
##                            x(k+1) = f(x(k)),
##
##   where f(x) is a continous function. Vector x is returned as output.
##
## Inputs
##   f     Function f(x)
##   N     Number of iterations
##   x1    Initial value of x_k
##   xmin  Minimum x in the cobwebbing plot
##   xmax  Maximum x in the cobwebbing plot
##
## Example
##
##   x=cobweb("max(1/4,x^2)",15,.1,0,pi);
##

## Author:   Simone Zuccher <zuccher@sci.univr.it>
## Created:  17 Apr 2008
## Modified: 


function x = cobweb(f,N,x1,xmin,xmax);

% Sanity check
if nargin < 5
   usage ("x = cobweb (f, N, x1, xmin, xmax)");
   return
endif

global func
func=f;

% Replace ^ * and / to make operations possible
func=strrep(func,"**","^");
func=strrep(func,"^",".^");
func=strrep(func,"*",".*");
func=strrep(func,"/","./");

% Number of iterations and initial condition
m=N;
x(1)=x1;

% Close all possible plots and change format to long exponential
close all
format long e

% Assign vx and vy needed for the cobwebbing plot
vx(1)=x(1);
vy(1)=0.0;

% Loop on all points
for n=1:1:m-1
  x(n+1)     = fk(x(n));
  vx(2*n)     = x(n);
  vy(2*n)     = x(n+1);
  vx(2*n+1)   = x(n+1);
  vy(2*n+1)   = x(n+1);
end

% Transpose to make it vector
x=x';

% Compute the arrows
m=length(vx);
for n=1:1:m-1
  arrx(n) = (vx(n)+vx(n+1))/2.-vx(n);
  arry(n) = (vy(n)+vy(n+1))/2.-vy(n);
end


% Plot x(k) versus k
figure(1);
__gnuplot_set__ auto
xlabel("k");ylabel("x_k");
plot(x,'*-r');legend('off');


% Second plot. Set legend back etc.
if exist("quiver") > 0
   figure(2);
endif
__gnuplot_set__ key bottom
leg1=["f(x) = " f];leg2="x";
xlabel("x_k");ylabel("");

% Vector t needed for final plot
t=linspace(xmin,xmax,500);


% If quiver exists, the cobwebbing plot is quite easy...
if exist("quiver") > 0
   quiver(vx(1:m-1),vy(1:m-1),arrx,arry);
   hold on
   plot(t,fk(t),[";" f ";"],t,t,";x;",vx,vy,'-r');
   hold off
   return
endif

% ...otherwise the cobwebbing plot is quite complicated...
M=[vx(1:m-1)' vy(1:m-1)' arrx' arry'];
ff=[t' (fk(t))'];
xx=[t' t'];
vv=[vx' vy'];
gp_cm=["__gnuplot_plot__ M t ''  w vector,\
                         vv t '' w l 1,\
	                 ff t 'f(x) = " f "' w l 2,\
			 xx t 'x' w l 3"];
eval(gp_cm);
return
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the function f(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fk=fk(x)
   global func
   eval(["fk=" func ";"]);
endfunction
