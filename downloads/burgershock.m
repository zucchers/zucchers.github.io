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

## Solve Riemann problem for Burger's equation
##
## u_t + u u_x = 0
##
## with IC
##
## u_l for x<0 and u_r for x>0.
## Comparisons are made between conservative and non conservative methods and 
## exact solution.
##

## Author:   Simone Zuccher <zuccher@sci.univr.it>
## Created:  10 Jun 2009
## Modified: 


clear all
close all
N=50;
xmin=-1;
xmax=1;
ul=1.2;
ur=0.4;

if(ul<ur)
   disp('Rarefaction')
else
   disp('Shock')
end

x=linspace(xmin,xmax,N);
u0=ur*ones(size(x));
u0(1:N/2)=ul;
um=(ul+ur)/2;

k=0.01;
tmax=0.6;
h=x(2)-x(1);
CFL=um*k/h;
if(abs(CFL)>1)
   error('CFL condition not satisfied');
end

unc = u0;
uc = u0;
i=0;
t=i*k;
while(t<tmax)
   unc = unc - k/h*unc.*(unc - [ul unc(1:end-1)]);
   uc  = uc  - k/h*(uc.^2/2 - [ul uc(1:end-1)].^2/2);
   i=i+1;
   t=i*k;
   plot(x,unc,'o-',x,uc,'*-',[xmin, x(2:end)+t*um],[ul, u0(2:end)],'-');
   axis([xmin xmax min(ul,ur)-.1*um max(ul,ur)+.1*um])
   title(['Equazione di Burger; h = ' num2str(h) ', '...
                               'k = ' num2str(k) ', '...
                               'CFL = ' num2str(CFL) ', '...
                               't = ' num2str(t)])
   legend('Non conservativo','Conservativo','Esatta','location','southwest')
   drawnow; 
   pause(0.2)
end
