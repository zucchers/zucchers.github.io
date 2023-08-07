% Author:   Simone Zuccher <zuccher@sci.univr.it>
% Created:  19 June 2010
% Last modified: 
%
% lbvp3x3 = Linear Boundary Value Problem of 3 equations in 3 unknowns
%
% Efficient script to solve on a generally uneven grid the linear ODE system
%
% v1'' + v2' + v1 + 3v3  = -sin y + 3 (sin y)^2 cos y
% v2'' - v2  - v1' + v3'= - 3 cos y + 2 sin y cos y - (sin y)^3
% v3'' + v2' = 2 (cos y)^3 - 7 cos y (sin y)^2 -sin y
%
% with BC
% @ y=0:    v1' = cos y,  v2  = cos y,  v3 = (sin y)^2 cos y
% @ y=ymax: v1  = sin y,  v2' = -sin y, v3 = (sin y)^2 cos y
%
% whose solution is (v1,v2,v3)=(sin y, cos y, (sin y)^2 cos y)
%
clear all
close all
%---------- Grid --------------------------------------------------------------%
N=500; ymax=2/3*pi; y=([0:(N-1)]./(N-1)).^2.5*ymax; %y=linspace(0,1,N)*ymax;
y=y(:)';
%---------- Constants ---------------------------------------------------------%
% equations and variables to make the code more readable
ne=3; e1=0; e2=1; e3=2; nv=3; v1=2; v2=0; v3=1;
% A temporary rectangular matrix which will become a sparse diagonal matrix
smin=-1; smax=1; sizeA=[N*ne,nv*(smax-smin+2)];
% md: column of rectangular matrix corresponding to the main diagonal
md = -(smin-1)*nv;
%---------- Variables ---------------------------------------------------------%
A=zeros(sizeA); rhs=zeros(sizeA(1),1);
%---------- Derivatives -------------------------------------------------------%
j=(2:N-1);
d1p(j)=1./(y(j+1)-y(j-1)); d1m(j)=-d1p(j);
d2p(j)=2*d1p(j)./(y(j+1)-y(j)); d2m(j)=2*d1p(j)./(y(j)-y(j-1));
d20(j)=-d2p(j)-d2m(j);
%---------- Matrix filling ----------------------------------------------------%
%
% First point (at y(1))
j=1; eq=(j-1)*ne+1;
A(eq+v1,md-e1+v1)=-1/(y(j+1)-y(j));
A(eq+v1,md-e1+v1+nv)=1/(y(j+1)-y(j));
A(eq+v2,md-e2+v2)=1;
A(eq+v3,md-e3+v3)=1;
rhs(eq+e1)=cos(y(j));
rhs(eq+e2)=cos(y(j));
rhs(eq+e3)=sin(y(j))^2*cos(y(j));
%
% Internal point from 2 to N-1
j=(2:N-1); eq=(j-1)*ne+1;
% First equation
A(eq+v1-nv,md-e1+v1-nv)=d2m(j);
A(eq+v2-nv,md-e1+v2-nv)=d1m(j);
A(eq+v1,md-e1+v1)=1+d20(j);
A(eq+v1+nv,md-e1+v1+nv)=d2p(j);
A(eq+v2+nv,md-e1+v2+nv)=d1p(j);
A(eq+v3,md-e1+v3)=3;
rhs(eq+e1)=-sin(y(j))+3*sin(y(j)).^2.*cos(y(j));
% Second equation
A(eq+v1-nv,md-e2+v1-nv)=-d1m(j);
A(eq+v2-nv,md-e2+v2-nv)=d2m(j);
A(eq+v3-nv,md-e2+v3-nv)=d1m(j);
A(eq+v2,md-e2+v2)=-1+d20(j);
A(eq+v1+nv,md-e2+v1+nv)=-d1p(j);
A(eq+v2+nv,md-e2+v2+nv)=d2p(j);
A(eq+v3+nv,md-e2+v3+nv)=d1p(j);
rhs(eq+e2)=-3*cos(y(j)) + 2*sin(y(j)).*(cos(y(j))).^2 - (sin(y(j))).^3;
% Third equation
A(eq+v2-nv,md-e3+v2-nv)=d1m(j);
A(eq+v3-nv,md-e3+v3-nv)=d2m(j);
A(eq+v3,md-e3+v3)=d20(j);
A(eq+v2+nv,md-e3+v2+nv)=d1p(j);
A(eq+v3+nv,md-e3+v3+nv)=d2p(j);
rhs(eq+e3)=2*(cos(y(j))).^3 -7*cos(y(j)).*(sin(y(j))).^2 -sin(y(j));
%
% Last point (at y(N))
j=N; eq=(j-1)*ne+1;
A(eq+v1,md-e1+v1)=1;
A(eq+v2,md-e2+v2)=1/(y(j)-y(j-1));
A(eq+v2-nv,md-e2+v2-nv)=-1/(y(j)-y(j-1));
A(eq+v3,md-e3+v3)=1;
rhs(eq+e1)=sin(y(j));
rhs(eq+e2)=-sin(y(j));
rhs(eq+e3)=sin(y(j))^2*cos(y(j));
%---------- Solution ----------------------------------------------------------%
A=spdiags(A,(1:sizeA(2))-md,sizeA(1),sizeA(1));
f=A\rhs;
var1=f(1+v1:nv:end);var2=f(1+v2:nv:end);var3=f(1+v3:nv:end);
%---------- Fancy plot --------------------------------------------------------%
k=1:5:N;
plot(y,sin(y),'-k',y(k),var1(k),'o',y,cos(y),'-k',y(k),var2(k),'*',...
     y,sin(y).^2.*cos(y),'-k',y(k),var3(k),'+')
legend('','v1','','v2','','v3','location','southwest')
xlabel('y');ylabel('v1, v2, v3')
title(['Linear boundary value problem for a system of 3 ODEs '...
       'with mixed boundary conditions'])
