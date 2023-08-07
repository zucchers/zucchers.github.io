clear all
close all
BLfp;
blasius;
close all

plot(UBL,yBL,'o',U,y,'-')
ylabel('y');
xlabel('u');
legend('BLfp','Blasius')

figure
plot(VBL,yBL,'o',V,y,'-')
ylabel('y');
xlabel('v');
legend('BLfp','Blasius')
