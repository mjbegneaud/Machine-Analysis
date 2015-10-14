%SCCA Cam Design
function [x y yprime ydblprime ytrplprime]=scca(method,riseorfall)
%
%Input variable rise or fall should be stated as 'rise' or 'fall'
%
%Input variable method should be given as one of the following:
%       'constant acceleration';
%       'modified trapezoid';
%       'simple harmonic';
%       'modified sine';
%       'cycloidal displacement';
%
switch lower(method)
    case('constant acceleration')
        b=0;
        d=0;
        c=1-b-d;
    case('modified trapezoid')
        b=0.25;
        d=0.25;
        c=1-b-d;
     case('simple harmonic')
        b=0;
        d=1;
        c=1-b-d;
     case('modified sine')
        b=0.25;
        d=0.75;
        c=1-b-d;
     case('cycloidal displacement')
        b=0.5;
        d=0.5;
        c=1-b-d;
end
%
Ca=4*pi^2/((pi^2-8)*(b^2-d^2)-2*pi*(pi-2)*b+pi^2);
Cv=Ca*((b+d)/pi+c/2);
Cj=Ca*pi/b;
%
x=0:0.01:1;
i=find(0<= x & x<=b/2);
x1=x(i);
i=find(b/2<x & x<=(1-d)/2);
x2=x(i);
i=find((1-d)/2<x & x<=(1+d)/2);
x3=x(i);
i=find((1+d)/2<x & x<=(1-b/2));
x4=x(i);
i=find((1-b/2)<x & x<=1);
x5=x(i);
clear i 
%Zone 1
y1=Ca*(b/pi*x1-(b/pi)^2*sin(pi/b*x1));
y1prime=Ca*(b/pi-b/pi*cos(pi/b*x1));
y1dblprime=Ca*sin(pi/b*x1);
y1trplprime=Ca*pi/b*cos(pi/b*x1);
%Zone 2
y2=Ca*(x2.^2/2+b*(1/pi-1/2)*x2+b^2*(1/8-1/pi^2));
y2prime=Ca*(x2+b*(1/pi-1/2));
y2dblprime=Ca*ones(size(x2));
y2trplprime= zeros(size(x2));
%Zone 3
y3=Ca*((b/pi+c/2)*x3+(d/pi)^2+b^2*(1/8-1/pi^2)-(1-d)^2/8-(d/pi)^2*cos(pi/d*(x3-(1-d)/2)));
y3prime=Ca*(b/pi+c/2+d/pi*sin(pi/d*(x3-(1-d)/2)));
y3dblprime=Ca*cos(pi/d*(x3-(1-d)/2));
y3trplprime=-Ca*pi/b*sin(pi/d*(x3-(1-d)/2));
%Zone 4
y4=Ca*(-x4.^2/2+(b/pi+1-b/2)*x4+(2*d^2-b^2)*(1/pi^2-1/8)-1/4);
y4prime=Ca*(-x4+b/pi+1-b/2);
y4dblprime=-Ca*ones(size(x4));;
y4trplprime=zeros(size(x4));
%Zone 5
y5=Ca*(b/pi*x5+2*(d^2-b^2)/pi^2+((1-b)^2-d^2)/4-(b/pi)^2*sin(pi/b*(x5-1)));
y5prime=Ca*(b/pi-b/pi*cos(pi/b*(x5-1)));
y5dblprime=Ca*(sin(pi/b*(x5-1)));
y5trplprime=Ca*pi/b*cos(pi/b*(x5-1));
x=[x1 x2 x3 x4 x5];
y=[y1 y2 y3 y4 y5];
yprime=[y1prime y2prime y3prime y4prime y5prime];
ydblprime=[y1dblprime y2dblprime y3dblprime y4dblprime y5dblprime];
ytrplprime=[y1trplprime y2trplprime y3trplprime y4trplprime y5trplprime];
switch lower(riseorfall)
    case ('fall')
        y=1-y;
        yprime=-yprime;
        ydblprime=-ydblprime;
        ytrplprime=-ytrplprime;
end
figure(1) 
subplot(2,2,1)
plot(x,y)
axis tight
xlabel('x')
ylabel('y')
title('Normalized S, V, A, & J During Rise Segment')
subplot(2,2,2)
plot(x,yprime)
axis tight
xlabel('x')
ylabel('y-prime')
subplot(2,2,3)
plot(x,ydblprime)
axis tight
xlabel('x')
ylabel('y-dble-prime')
subplot(2,2,4)
plot(x,ytrplprime)
axis tight
xlabel('x')
ylabel('y-trple-prime')











