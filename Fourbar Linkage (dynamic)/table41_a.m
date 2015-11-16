% Script used to calculate the dynamic forces present in a fourbar linkage throughout its range of motion.
% Author: Matthew J. Begneaud
% In collaboration with:
%   Ronald Kisor
%   Chandler Lagarde
%   Jace Delcambre
% Date: 11/9/15

% This routine solves for the unknown thetas, theta3, theta4, omega3,
% omega4, alpha3, and alpha 4 for a fourbar linkage.
% The dynamic forces and torque are then solves for using a mass density of
% 1kg/25cm
clear
clc

% This allows user to set the values of these variables in functions that 
% are called without having to pass them as input parameters
global a b c d t2 w2 a2

% Angular range of motion (in deg) for link2 (Theta2 range)
t2start=0; 
t2stop=360; 

% Initial guess (in deg) for Theta3 and Theta4 (vector of guesses included)
t3est=106.6;
t4est=131.8;
tunknownguess=[t3est; t4est]*pi/180; 

% Link lengths in cm
a=2;
b=7;
c=9;
d=6;

% Omega2 and Alpha2
w2=10;
a2=0;


% Setup for step size of Theta2
t2=t2start*pi/180;
n=floor(((t2stop+5)-t2start)/5); 

% Start counter
i=1;


% %%%%%%%%%%%%%%% Kinematic Analysis %%%%%%%%%%%%%%%%%%%%%

% Allocate memory space for unknownthetas matrix and Fval matrix
unknownthetas=zeros(n,2); 
Fval=zeros(n,2); 


% Main loop for Thetas
while (1)  
    [unknownthetas(i,:),Fval(i,:)]=fsolve(@thetwoeqtns,tunknownguess);
    t2=t2+5*pi/180;
    tunknownguess=unknownthetas(i,:);
    i=i+1;
    if i>n
        break
    end
end
% Generate full vector of Theta2, Theta3, and Theta4 values
t2=(t2start:5:t2stop)'*pi/180; 
t3 = unknownthetas(:,1);
t4 = unknownthetas(:,2);


% Main loop for Omegas
omegas = zeros(2,n);
for i = 1:n
    A = [-b*sin(t3(i)), c*sin(t4(i));...
        b*cos(t3(i)), -c*cos(t4(i))];
    rhs = a*w2*[sin(t2(i)), -cos(t2(i))]';
    omegas(:,i) = A\rhs;
end
% Generate full vector of Omega3 and Omega4 values
w3 = omegas(1,:)';
w4 = omegas(2,:)';


% Main loop for Alphas
alphas = zeros(2,n);
for i = 1:n
    B = [-b*sin(t3(i)), c*sin(t4(i));...
        b*cos(t3(i)), -c*cos(t4(i))];
    arhs = a*a2*[sin(t2(i)), -cos(t2(i))]' + a*(w2^2)*[cos(t2(i)), sin(t2(i))]'...
        + b*(w3(i)^2)*[cos(t3(i)), sin(t3(i))]' - c*(w4(i)^2)*[cos(t4(i)), sin(t4(i))]';
    alphas(:, i) = B\arhs;
end
% Generate full vector of Alpha3 and Alpha4 values
a3 = alphas(1,:)';
a4 = alphas(2,:)';


% % Plot theta3 versus theta2 and theta4 versus theta2 in separate plots
% subplot(2,1,1)
% plot(t2*180/pi,unknownthetas(:,1)*180/pi)
% axis tight
% grid on
% xlabel('Theta2, deg')
% ylabel('Theta3, deg')
% subplot(2,1,2)
% plot(t2*180/pi,unknownthetas(:,2)*180/pi)
% axis tight
% grid on
% xlabel('Theta2, deg')
% ylabel('Theta4,deg')

% % Plot omega3 vs theta2 | omega4 vs theta2 | alpha3 vs theta2 |
% % alpha4 vs theta2
% subplot(2, 2, 1)
% plot(t2*180/pi, w3)
% xlabel('Theta2, deg')
% ylabel('Omega3, rad/s')
% axis tight
% grid on
% subplot(2,2,2)
% plot(t2*180/pi, w4)
% xlabel('Theta2, deg')
% ylabel('Omega4, rad/s')
% axis tight
% grid on 
% subplot(2,2,3)
% plot(t2*180/pi, a3)
% xlabel('Theta2, deg')
% ylabel('Alpha3, rad/s^2')
% axis tight
% grid on
% subplot(2,2,4)
% plot(t2*180/pi, a4)
% xlabel('Theta2, deg')
% ylabel('Alpha4, rad/s^2')
% axis tight
% grid on



% %%%%%%%%%%%%%%%% Dynamic Analysis %%%%%%%%%%%%%%%%%%


% Center of gravity accelerations & inertia values
% link 2
m2 = (a/25);
i2 = (1/3)*m2*(a^2);
axg2 = -(a/2)*w2^2*cos(t2)-(a/2)*a2*sin(t2);
ayg2 = -(a/2)*w2^2*sin(t2)+(a/2)*a2*cos(t2);
% link 3
m3 = (b/25);
i3 = (1/12)*m3*(b^2);
axg3 = a*-cos(t2)*w2^2+a*-sin(t2)*a2 + (b/2)*-cos(t3).*w3.^2+(b/2)*-sin(t3).*a3;
ayg3 = a*-sin(t2)*w2^2+a*cos(t2)*a2 + (b/2)*-sin(t3).*w3.^2+(b/2)*cos(t3).*a3;
% link 4
m4 = (c/25);
i4 = (1/3)*m4*(c^2);
axg4 = (c/2)*-cos(t4).*w4.^2+(c/2)*-sin(t4).*a4;
ayg4 = (c/2)*-sin(t4).*w4.^2+(c/2)*cos(t4).*a4;



% Main loop for forces and torques
for i=1:n
    dynamicmatrix = [1,0,0,a*sin(t2(i)),a*cos(t2(i)),0,0,0,0;...
    0,1,0,-1,0,0,0,0,0;...
    0,0,1,0,1,0,0,0,0;...
    0,0,0,1,0,-1,0,0,0;...
    0,0,0,0,-1,0,1,0,0;...
    0,0,0,(b/2)*sin(t3(i)),(b/2)*cos(t3(i)),(b/2)*sin(t3(i)),(b/2)*cos(t3(i)),0,0;...
    0,0,0,0,0,-c*cos(t4(i)),b*sin(t4(i)),0,0;...
    0,0,0,0,0,1,0,-1,0;...
    0,0,0,0,0,0,-1,0,1];
    dynamicrhs = [i2*a2;m2*axg2(i);m2*ayg2(i);m3*axg3(i);m3*ayg3(i);i3*a3(i);i4*a4(i);m4*axg4(i);m4*ayg4(i)];   
    unknown_forces(:,i) = dynamicmatrix\dynamicrhs;
end
% Torque and force values
torque = unknown_forces(1,:);
F1x = unknown_forces(2,:);
F1y = unknown_forces(3,:);
F2x = unknown_forces(4,:);
F2y = unknown_forces(5,:);
F3x = unknown_forces(6,:);
F3y = unknown_forces(7,:);
F4x = unknown_forces(8,:);
F4y = unknown_forces(9,:);


