clear all
close all
% 0 < b < 1 : origin stable
disp('For 0<b<b*=1 the origin is the only equilibrium and is stable...')
[t,x]= mylorenz(0.01,5,[10 1/2 8/3],[1 15 1],0);
disp('Press enter to continue...'); pause;
[t,x]= mylorenz(0.01,5,[10 1/2 8/3],[-15 -5 8],1);
disp('Press enter to continue...'); pause;
[t,x]= mylorenz(0.01,5,[10 1/2 8/3],[15 15 15],1);
% 1 < b < b** : origin unstable + 2 more stable equilibria
disp('For b*<b<b** there are 3 equilibria, the 2 new are stable...')
disp('Press enter to continue...'); pause;
[t,x]= mylorenz(0.05,10,[10 3/2 8/3],[-0.3 -0.7 0.5],0);
%disp('Press enter to continue...'); pause;
[t,x]= mylorenz(0.05,10,[10 3/2 8/3],[-0.2 0.3 -0.5],1);
%disp('Press enter to continue...'); pause;
[t,x]= mylorenz(0.05,10,[10 3/2 8/3],[.01 0.01 .01],1);
[t,x]= mylorenz(0.05,10,[10 3/2 8/3],[-.01 -0.01 -.01],1);
[t,x]= mylorenz(0.05,10,[10 3/2 8/3],[1 -.5 -.5],1);
[t,x]= mylorenz(0.05,10,[10 3/2 8/3],[1 .5 .5],1);
[t,x]= mylorenz(0.05,10,[10 3/2 8/3],[0 0 1],1);
[t,x]= mylorenz(0.05,10,[10 3/2 8/3],[0 0 -1],1);
disp('For b ~= b** the 2 new equilibria are stable but it takes a while...')
disp('Press enter to continue...'); pause;
[t,x]= mylorenz(0.01,110,[10 22.5 8/3],[-5 20 22],0);
disp('Press enter to continue...'); pause;
figure(3)
plot(t,x(:,1))
disp('For b>b** the 2 new equilibria are unstable...')
disp('Press enter to continue...'); pause;
x0=[3 15 1];
[t,x]= mylorenz(0.01,50,[10 28 8/3],x0,0);
disp('For b>b** check out the sensitivity to initial conditions...')
disp('Press enter to continue...'); pause;
[t,xp]= mylorenz(0.01,50,[10 28 8/3],1.00001*x0,0);
%disp('Press enter to continue...'); pause;
figure(3)
plot(t,x(:,1),t,xp(:,1))
disp('For b=100.55 there is one stable periodic orbit...')
disp('Press enter to continue...'); pause;
[t,x]= mylorenz(0.005,20,[10 100.55 8/3],[-1.64 -4.9 64.13],0);
figure(3)
plot(x(:,1),x(:,3))
disp('For b=99.8 you can see the period doubbling...')
disp('Press enter to continue...'); pause;
[t,x]= mylorenz(0.005,20,[10 99.8 8/3],[-1.64 -4.9 64.13],1);
figure(3)
plot(x(:,1),x(:,3))
