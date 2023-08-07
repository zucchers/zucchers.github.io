function [t, X] = planarsys(f,g,xmin,xmax,ymin,ymax,dt,tmax,x0,y0,rep);
% This version has been tested on:
% + GNU Octave 3.0
% + Matlab 7.4.0.336 (R2007a)
%
% planarsys(f,g,xmin,xmax,ymin,ymax,dt,tmax,x0,y0,rep)
%
%   Given the ODE system
%      _
%     |
%     |  d x
%     |  --- = f(x,y)
%     |  d t
%    <
%     |  d y
%     |  --- = g(x,y)
%     |  d t
%     |_        
%
%
%   planarsys plots:
%      (a) the vector field if quiver exists 
%      (b) the horizontal and vertical isoclines 
%      (c) the equilibria (and show them on the screen)
%      (d) the trajectory x(t),y(t) obtained from the initial condition (x0,y0)
%          either marching forward or backward in time (see dt and tmax)
%
%
% Inputs
%   f     Function f(x,y)
%   g     Function g(x,y)
%   xmin  Minimum x in the phase portrait plot
%   xmax  Maximum x in the phase portrait plot
%   ymin  Minimum y in the phase portrait plot
%   ymax  Maximum y in the phase portrait plot
%   dt    Time-step for time integration (can be negative for backward-in-time
%         integration)
%   tmax  Maximum time (must be negative if dt < 0 )
%   x0    Initial condition x(0)=x0
%   y0    Initial condition y(0)=y0
%   rep   Flag: 1 - replot the new trajectory on the existing plot
%               0 - plot everything from scratch
%
%
% Outputs
%   t     Vector of time
%   X     Solution, i.e. x=X(:,1) and y=X(:,2)
%
% Examples
%
% 1. Plot the vector field, equilibria and a trajectory from scratch:
%
%    [t, X] = planarsys('x*(5-4*x-y)','y*(4-2*x-3*y)',-.2,1.6,-.2,1.6,.05,4,.01,.1,0);
%
% 2. Plot a new trajectory on the existing plot:
%
%    [t, X] = planarsys('x*(5-4*x-y)','y*(4-2*x-3*y)',-.2,1.6,-.2,1.6,.05,4,.01,1.5,1);
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
% Created:  14 Jun 2008


% Sanity check
if nargin < 11
   usage ('planarsys(f,g,xmin,xmax,ymin,ymax,dt,tmax,x0,y0,rep)');
   return
end

% Check if quiver exists...
if((exist('quiver')==0) )
   disp('WARNING: quiver was not found, therefore vector field is not shown.');
   ifquiver = false;
else
   ifquiver=true;
end

% Check whether Octave or Matlab is running...
if(exist('octave_config_info')>0)
   ifOctave = true;
   ifMatlab = false;
else
   ifOctave = false;
   ifMatlab = true;
end


%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%
%%% Options changeable by the user
%%%
usequiver=true; % if TRUE use quiver if available, otherwise don't
axn   = 25;     % Number of arrows in x (vector field)
ayn   = 25;     % Number of arrows in y (vector field)
anorm = true;   % if TRUE vectors are normalized (=all same length) 
                % if FALSE vectors have different lenghts
xn    = 100;    % Number of points in x for plotting f(x,y)=0 and g(x,y)=0
yn    = 100;    % Number of points in y for plotting f(x,y)=0 and g(x,y)=0
% Equilibrium point(s) need guees to be found:
ngord = 20;     % Number of guesses when |f(x,y)| + |g(x,y)| is ordered
ng    = 5;      % Number of guesses in x and y when scanning plot area (ng x ng)
                % NOTE: the total number of guesses is ngord + ng^2
%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%






% This is needed for the inline function F
f=strrep(f,'x','x(1)');
f=strrep(f,'y','x(2)');
g=strrep(g,'x','x(1)');
g=strrep(g,'y','x(2)');
% This is for Octave
fstro = [ '[ ' f '  ' g  ']' ]; 
Fo = inline(fstro, 'x');
% This is for Matlab
fstrm = [ '[ ' f ' ; ' g  ']' ]; 
Fm = inline(fstrm, 'x');

% If replot flag is true, close all windows and start over from scratch   
if(rep ~= 1)
   % Close all open plots
   close all
   [x,y]=meshgrid(linspace(xmin,xmax,xn),linspace(ymin,ymax,yn));
   dx=zeros(size(x));
   dy=zeros(size(x));
   for i=1:xn
      for j=1:yn
         tmp=Fm([x(i,j);y(i,j)]);
	 dx(i,j)=tmp(1);
	 dy(i,j)=tmp(2);
      end
   end
   % Plot horizontal and vertical isoclines
   contour(x,y,dy,[0 0],'-r')
   hold on
   contour(x,y,dx,[0 0],'-g')
   hold on
   axis([xmin xmax ymin ymax]);


   % Look for equilibrium point(s)
   % Strategy: organize a matrix aa = |f(x,y)| + |g(x,y)|
   % and extract the first ngord points (x,y) ordered from the smallest 
   % element of aa
   aa=abs(dx)+abs(dy);            % aa = |f(x,y)| + |g(x,y)|
   [aord iaord] = sort(aa(:));    % sort everything out
   xg = x(iaord(1:ngord));
   yg = y(iaord(1:ngord));

   % Compute equilibria numerically
   [xgt,ygt]=meshgrid(linspace(xmin,xmax,ng),linspace(ymin,ymax,ng));
   
   xg=[xg; xgt(:)];
   yg=[yg; ygt(:)];
%  New size of guess points
   ng=length(xg);
   for i=1:ng
      % Solve the nonlinear system (different options in Matlab)
      if(ifOctave)
         [eqt,info,exf] = fsolve (Fm, [xg(i); yg(i)]);
      else
	 [eqt,info,exf] = fsolve (Fm, [xg(i); yg(i)],...
	                        optimset('Display','off'));
      end
      if(i==1)
         eq=eqt';
      else   
         eq = [eq;eqt'];
      end   
   end 

   % Make equilibria complex to sort them easier
   leq=length(eq);
   eqc=eq(:,1)+sqrt(-1)*eq(:,2);
   eqc=str2num(num2str(eqc,'%16.4f'));
   eqc=sort(eqc);
   % This is to get rid of too much precision
   % Get back to real x and y
   eq(:,1)=real(eqc);
   eq(:,2)=imag(eqc);
   % Find out the non-repeated equilibria 
   b=eq;
   a = diff(b);
   a = str2num(num2str(a,'%16.4f'));
   c=abs(a(:,1))+abs(a(:,2));
   [a i] = find(c);
   a = [0 a']'+1;
   % Set eq to the correct values sorted by their distance from the origin
   eq=b(a,:);
   % Show equilibria on the screen
   seq=size(eq);
   disp(['Found ' num2str(seq(1)) ' equilibrium point(s):']);
   format % To keep only 5 significant digits
   disp(eq)
   disp('If you forsee more equilibria, please change ngord or ng in planarsys.m')
   % ...and finally plot them 
   hold on
   plot(eq(:,1),eq(:,2),'or');
   % Done with equilibrium points
   
   % Vector field (only if quiver exists) 
   if(ifquiver && usequiver) 
      [x,y]=meshgrid(linspace(xmin,xmax,axn),linspace(ymin,ymax,ayn));
      dx=zeros(size(x));
      dy=zeros(size(x));
      for i=1:axn
	 for j=1:ayn
	    tmp=Fm([x(i,j);y(i,j)]);
	    dx(i,j)=tmp(1);
	    dy(i,j)=tmp(2);
	 end
      end
      S=.5;
      if(ifOctave)
	 S=1;
      end
      L = min(xmax-xmin,ymax-ymin)/3*sqrt(dx.^2 + dy.^2);
      % Set L=1 for unscaled vector field
      if(~anorm)
	 L=1;
      end   
      % Plot the flow field and the equilibria if quiver exists... 
      hold on
      quiver(x,y,dx./L,dy./L,S)
      hold on
      axis([xmin xmax ymin ymax]);
   end % End of quiver part
end % End of replot part

% Re-define F as a function of t and x and solve the ODE system 
Fm = inline(fstrm, 't', 'x');
if (exist('lsode')>0)
% Octave
   t = 0:dt:tmax;
   X =  lsode(Fo, [x0 y0], t);
elseif(exist('ode45')>0)
% Matlab
   [t,X]=ode45(Fm,[0,tmax],[x0;y0]);
end

% Plot intial condition
hold on
plot(X(1,1), X(1,2), 'ok')
% Plot solution
hold on
plot(X(:,1), X(:,2), '-k')
% Set labels and stuff
xlabel('x(t)');
ylabel('y(t)');
axis([xmin xmax ymin ymax]);
larrow=min(xmax-xmin,ymax-ymin)/20;
if(dt>0)
   iarrow=10;
else
   iarrow=length(t)-10;  
end
if(iarrow>length(t))
   iarrow = floor(length(t)/2);
end
myarror(t,X,iarrow,larrow,20)
return
% End of main planarsys


function myarror(t,X,i,l,alpha)
alpha=alpha/180*pi;
signdt=sign(t(i)-t(i-1));
l=l*signdt;
beta=atan2(X(i,2)-X(i-1,2),X(i,1)-X(i-1,1));
hold on
plot([X(i,1) X(i,1)-l*cos(alpha+beta)],[X(i,2) X(i,2)-l*sin(alpha+beta)],'-k')
hold on
plot([X(i,1) X(i,1)-l*cos(beta-alpha)],[X(i,2) X(i,2)-l*sin(beta-alpha)],'-k')
