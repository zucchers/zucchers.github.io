%% Copyright (C) 2008 Simone Zuccher
%% Purpose:     To show how gps works
%% Author:      Simone Zuccher
%% Created:     01 Nov 2008
%% Suggestions: zuccher@sci.univr.it

% Make sure previous gps windows is closed
gps

% Close all previous stuff
clear all
close all

% Clear the terminal window
clc

% Delay between plots
delay=3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Simple gnuplot commands
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A first simple gnuplot command
gps("plot sin(x) t '$\\sin x$' w lp pt 7;\
     set auto;\
     set yrange[-1:1.6];\
     set xlabel '$x$';\
     set title 'Let''s plot a sine with bullets (filled circles)'")
disp(["wait please " num2str(delay) " seconds..."]);fflush(stdout);pause(delay);
% Add a cosine
gps("replot cos(x) t '$\\cos x$' w lp pt 4;\
     set title 'Let''s add a cosine (empty squares) to the sine (filled circles)'")
% This makes the png figure you see on the webpage
gps("set terminal png transparent;set output 'demogps2d1.png';rep;set term x11")
% This makes the ps figure you see on the webpage
gps("ps","demogps2d1")
disp(["wait please " num2str(delay) " seconds..."]);fflush(stdout);pause(delay);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Stable spiral
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close previous window
gps
% Window ranges
xmin=-2.0;
xmax=2.;
ymin=-2.0;
ymax=2.;

% Model constant
global aa;

% Set matrix
aa=[-1 1 -2 -1];

% Set initial conditions
x0=[-1 2; 0 2; 1 2; 1.7 2; 2 1; 2 0; 2 -1];
x0= [x0; [-x0(:,1) -x0(:,2)]];

% Set time options and create the vector
tmax=5;
tmin=0;
dt=.005;
t = tmin:dt:tmax;

% Definition of the dynamical system
function xdot=dsys(x, t)
  global aa;
  u = x(1);
  v = x(2);
  xdot(1) = aa(1)*u + aa(2)*v;
  xdot(2) = aa(3)*u + aa(4)*v;
endfunction

% Prepare the plot area
gps(["set size square;\
      set xrange[" num2str(xmin) ":" num2str(xmax) "];\
      set yrange[" num2str(ymin) ":" num2str(ymax) "];\
      unset key;\
      set xlabel '$x(t)$';\
      set ylabel '$y(t)$';\
      set title 'This is a stable spiral';\
      " ])

% This is to generate the spirals
for i=1:length(x0)
   % Solve the ode
   x =  lsode("dsys", x0(i,:), t)';
   if (i>1)
      pltoption="replot";
   else
      pltoption="plot";
   endif
   % Plot the current spiral
   gps(pltoption,[x(1,:)' x(2,:)'],"t '' w l 1")
   % This is for the arrow
   it=85;
   scale=0.05*max(abs(xmax-xmin),abs(ymax-ymin));
   % Plot the arrow for the current spiral
   gps(["set arrow from " num2str(x(1,it)) "," num2str(x(2,it)) " to "\
             num2str(x(1,it+1)) "," num2str(x(2,it+1)) " size " num2str(scale)\
	     ", 20 lt 1"]);
end
% Add a filled circle in the origin
gps("replot",[0 0],"t '' w p pt 7 ps 2")
% Generate vecotr field
yn=xn=20;
[X,Y]=meshgrid(linspace(xmin,xmax,xn),linspace(ymin,ymax,yn));
dX=aa(1)*X + aa(2)*Y; 
dY=aa(3)*X + aa(4)*Y;
% Normalize vectors. 
L = 7*sqrt(dX.^2 + dY.^2);
% To avoid normalization leave the next line uncommented, otherwise comment it
L = .5*max(max(L))*ones(size(dX));
% Finally, plot. Note that this single line is equivalent to quiver
gps("replot",[X(:) Y(:) (dX./L)(:) (dY./L)(:)],"w vectors lt 1")
% Remove the title just for the plots
gps("set title ''")
% This makes the ps figure you see on the webpage
gps("psbk","demogps2d2")
% This makes the png figure you see on the webpage
gps("set terminal png transparent;set output 'demogps2d2.png';rep;set term x11")
% Set the title back
gps("set title 'This is a stable spiral'");
disp(["wait please " num2str(delay) " seconds..."]);fflush(stdout);pause(delay);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Sombrero
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close previous window
gps
% Number of points in x and y
n=41;
% Prepare data for plot
tx = ty = linspace (-8, 8, n)';
[xx, yy] = meshgrid (tx, ty);
r = sqrt (xx .^ 2 + yy .^ 2) + eps;
tz = sin (r) ./ r;
% Plot using surface and various gnuplot stuff
gps("splots",tx,ty,tz,"w l t '';\
        set surface;\
        set hidden3d;\
        set ticslevel 0;\
        set xlabel '$x$';\
        set ylabel '$y$';\
        set label '$z=\\sin\\left(\\sqrt{x^2+y^2}\\right)$' at  0,0,1.1;\
        set title 'This is the famous sombrero';\
        ")
% Remove the title just for the plots
gps("set title ''")
% This makes the ps figure you see on the webpage
gps("psbk3d","demogps3d1")
% This makes the png figure you see on the webpage
gps("set terminal png transparent;set output 'demogps3d1.png';rep;set term x11")
gps("set title 'This is the famous sombrero'");
disp(["wait please " num2str(delay) " seconds..."]);fflush(stdout);pause(delay);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Lorenz
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close previous window
gps

% Model constant
global aa;

% Set model parameters
aa=[10 28 8./3.];

% Equilibrium points
eq = [ 0 0 0];
eq = [eq; sqrt(aa(2)*aa(3)-aa(3)) sqrt(aa(2)-1)*sqrt(aa(3)) aa(2)-1];
eq = [eq; -sqrt(aa(2)*aa(3)-aa(3)) -sqrt(aa(2)-1)*sqrt(aa(3)) aa(2)-1];

% Set initial condition
x0=[3 15 1];

% Set time parameters
tmax=50;
tmin=0;
dt=.01;
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

i=1;
% Solve the ode
x =  lsode("dsys", x0(i,:), t');
% Plot the solution as a single line in 3D
gps("splot",[x(:,1) x(:,2) x(:,3)],"t '' w l 1")
% Plot an arrow
it=10;
gps(["set arrow from " num2str(x(it,1)) "," num2str(x(it,2)) \
     "," num2str(x(it,3)) " to "\
    num2str(x(it+1,1)) "," num2str(x(it+1,2)) "," \
    num2str(x(it+1,3))  " size graph 0.02, 20 lt 1"]);
% Plot another one
it=3030;
gps(["set arrow from " num2str(x(it,1)) "," num2str(x(it,2)) \
     "," num2str(x(it,3)) " to "\
    num2str(x(it+1,1)) "," num2str(x(it+1,2)) "," \
    num2str(x(it+1,3))  " size graph 0.02, 20 lt 1"]);
% Set the axis stuff
gps(["set auto;\
      unset key;\
      set xlabel '$x(t)$';\
      set ylabel '$y(t)$';\
      set zlabel '$z(t)$' offset graph 0.05,0.05,.55;\
      set ticslevel 0;\
      set title 'This is the famous Lorenz attractor';\
      " ])
% Plot the stable and unstable equilibria
gps("replot",[0 0 0],"t '' w p pt 6 ps 2",\
            [sqrt(aa(2)*aa(3)-aa(3)) sqrt(aa(2)-1)*sqrt(aa(3)) aa(2)-1],\
	           "t '' w p pt 7 ps 2",\
            [-sqrt(aa(2)*aa(3)-aa(3)) -sqrt(aa(2)-1)*sqrt(aa(3)) aa(2)-1],\
	           "t '' w p pt 7 ps 2")
% Remove the title just for the plots
gps("set title ''")
% This makes the png figure you see on the webpage
gps("psbk3d","demogps3d2")
% This makes the png figure you see on the webpage
gps("set terminal png transparent;set output 'demogps3d2.png';rep;set term x11")
% Set back the title
gps("set title 'This is the famous Lorenz attractor'")
