% Matthew J. Begneaud, 10/12/15
% Run SCCA design code
% apply the scale factors h, omega, and beta
%
% This is a rise dwell fall dwell (RDFD) cam
% The rate of rotation is 10 RPM

clear
clc

% %%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%

h=1.5; %units are mm
omega=2*pi*10/60; 
betarise=45*pi/180;
highdwell=150*pi/180;
betafall=90*pi/180;
lowdwell=2*pi-(betarise+highdwell+betafall);
%


N=100; 
        
        

% %%%%%%%%%%%%%%%%%%%% SCCA calls %%%%%%%%%%%%%%%%%%%%%


[xrise yrise yprise ydblprise ytrplprise]=scca('modified sine','rise');




% %%%%%%%%%%%%%%%%%% Reassignment for Fall %%%%%%%%%%%%%%%%%%%%%%%
xfall=xrise;
yfall=1-yrise;
ypfall=-yprise;
ydblpfall=-ydblprise;
ytrplpfall=-ytrplprise;




% %%%%%%%%%%%%%%%%%% Y values for Dwells %%%%%%%%%%%%%%%%%%%%%%%%%%%

yhighdwell=ones(1,N);
yphighdwell=zeros(1,N);
ydblphighdwell=zeros(1,N);
ytrplphighdwell=zeros(1,N);

ylowdwell=zeros(1,N);
yplowdwell=zeros(1,N);
ydblplowdwell=zeros(1,N);
ytrplplowdwell=zeros(1,N);


% %%%%%%%%%%%%%%%%%%%% Theta Values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The theta vectors corresponding to the rise, highdwell, fall, and low
%dwell segments are constructed below
thetarise=xrise*betarise;
thetahighdwell=betarise+[1:N]*highdwell/N;
thetafall=max(thetahighdwell)+betafall*xfall;
thetalowdwell=max(thetafall)+[1:N]*lowdwell/N;


% %%%%%%%%%%%%%%%%%%%%%%% Generated Vectors %%%%%%%%%%%%%%%%%%%%
% Below, the entire 360 degree set of theta, S, V, A, and J values are
% assembled into vectors
theta=[thetarise thetahighdwell thetafall thetalowdwell];
S=[yrise yhighdwell yfall ylowdwell]*h;
V=omega*h*[yprise/betarise yphighdwell ypfall/betafall yplowdwell];
A=omega^2*h*[ydblprise/betarise^2 ydblphighdwell ydblpfall/betafall^2 ydblplowdwell];
J=omega^3*h*[ytrplprise/betarise^3 ytrplphighdwell ytrplpfall/betafall^3 ytrplplowdwell];



% %%%%%%%%%%%%%%%%%%%% Plot Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 1 is generated by the SCCA function and shows only the normalized
%lift, velocity, and acceleration during the rise segment for the selected
%SCCA function type

figure(2)
subplot(2,2,1)
plot(theta*180/pi,S)
axis tight
xlabel('Theta, deg')
ylabel('Displacement, mm')
grid on
subplot(2,2,2)
plot(theta*180/pi,V)
axis tight
xlabel('Theta, deg')
ylabel('Velocity, mm/s')
grid on
subplot(2,2,3)
plot(theta*180/pi,A)
xlabel('Theta, deg')
ylabel('Acceleration, mm/s^2')
axis tight
grid on
subplot(2,2,4)
plot(theta*180/pi,J)
axis tight
xlabel('Theta, deg')
ylabel('Jerk, mm/s^3')
grid on
maxS=max(S)
maxV=max(V) 
maxA=max(A) 
maxJ=max(J)

