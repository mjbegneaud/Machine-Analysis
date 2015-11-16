% Script used for kinematic and dynamic analysis of an offset slider mechanism.
% Author: Matthew Begneaud
% In collaboration with:
%     Jace Delcambre  
%     Chris Falcon
%     Chandler Lagarde
% Date: 11/16/15

clear 
clc

% Angular range of motion (in deg) for link2 (Theta2 range)
t2start=0; 
t2stop=360; 
t2=t2start*pi/180;
n=floor(((t2stop+5)-t2start)/5); 
t2=(t2start:5:t2stop)*pi/180; 

% Link Lengths (m)
a = 3*0.01;
b = 13*0.01;
c = 0*0.01;
% Radii of Gyration (cm)
rg2 = 1*0.01;
rg3 = 6.5*0.01;

% Mass of links (kg)
m2 = 2.5;
m3 = 1 ;
mb = 0.5;

% Coeff. of Friction
u = 0.1;

% Omega2 and Alpha2
w2 = 2000*2*pi/60;
a2 = 0;



% %%%%%%%%%%%%%%% Kinematic Analysis %%%%%%%%%%%%%%%%

% Position
t3 = zeros(1,n);
for i=1:1:n
    t3(:,i) = asin((a*sin(t2(i))-c)/b);
end

% Velocity
velocities = zeros(2,n);
for i=1:1:n
    V = [1,-b*sin(t3(i)) ; 0, b*cos(t3(i))];
    Vrhs = a*w2*[-sin(t2(i));cos(t2(i))];
    velocities(:,i) = V\Vrhs;
end
ddot = velocities(1,:);
w3 = velocities(2,:);

% Acceleration
accelerations = zeros(2,n);
for i=1:1:n
   A = [1,-b*sin(t3(i)) ; 0, b*cos(t3(i))];
   Arhs = a*a2*[-sin(t2(i));cos(t2(i))]-a*w2.^2*[cos(t2(i));sin(t2(i))] + b*w3(i).^2*[cos(t3(i));sin(t3(i))];
   accelerations(:,i) = A\Arhs;
end
ddbldot = accelerations(1,:);
a3 = accelerations(2,:);



% %%%%%%%%%%%%%%%%% Dynamic Analysis %%%%%%%%%%%%%%%%%%%

% Center of Gravity Accelerations
% link 2
xdbldg2 = -rg2*(a2*sin(t2) + w2^2*cos(t2));
ydbldg2 = rg2*(a2*cos(t2) - w2^2*sin(t2));
i2 = m2*rg2^2+m2*rg2^2; % I believe this is incorrect, but it gets multiplied by 2 when used, so not important in this case...

% link 3
xdbldg3 = ddbldot - a3*(b-rg3).*sin(t3) - w3.^2*(b-rg3).*cos(t3);
ydbldg3 = 0 + a3*(b-rg3).*cos(t3) - w3.^2*(b-rg3).*sin(t3);
i3 = m3*rg3^2;


% Force calculations
% Sign of ddot is required for determining direction of friciton force on
% slider.
sign_ddot = zeros(1,n);
for i =1:1:n
    sign_ddot(:,i) = ddot(i)/abs(ddot(i));
end
sign_ddot(1)=0;

% Calculation of dynamic equations using matrices
for i=1:1:n
   dynamicmatrix = [1,0,-1,0,0,0,0,0;...
       0,1,0,-1,0,0,0,0;...
       0,0,a*sin(t2(i)),-a*cos(t2(i)),0,0,0,1;...
       0,0,1,0,-1,0,0,0;...
       0,0,0,1,0,1,0,0;...
       0,0,-rg3*sin(180-t3(i)),-rg3*cos(180-t3(i)),-(b-rg3)*sin(180-t3(i)),(b-rg3)*cos(180-t3(i)),0,0;...
       0,0,0,0,1,0,-u*sign_ddot(i),0;...
       0,0,0,0,0,-1,1,0];
   dynamicrhs = [m2*xdbldg2(i);m2*ydbldg2(i);i2*a2;m3*xdbldg3(i);m3*ydbldg3(i);i3*a3(i);mb*ddbldot(i);0];
   unknown_forces(:,i) = dynamicmatrix\dynamicrhs;
end
% Torque and force values
tau = unknown_forces(8,:);
F1x = unknown_forces(1,:);
F1y = unknown_forces(2,:);
F2x = unknown_forces(3,:);
F2y = unknown_forces(4,:);
F3x = unknown_forces(5,:);
F3y = unknown_forces(6,:);
N = unknown_forces(7,:);
% Resultant forces
F1 = sqrt(F1x.^2+F1y.^2);
F2 = sqrt(F2x.^2+F2y.^2);
F3 = sqrt(F3x.^2+F3y.^2);



% %%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%
% Plot tau
subplot(2,2,1)
plot(t2*180/pi, tau)
xlabel('Theta2, deg')
ylabel('Tau, N-m')
axis tight
grid on
% Plot F1
subplot(2,2,2)
plot(t2*180/pi, F1)
xlabel('Theta2, deg')
ylabel('F1, N-m')
axis tight
grid on
% Plot F2
subplot(2,2,3)
plot(t2*180/pi, F2)
xlabel('Theta2, deg')
ylabel('F2, N-m')
axis tight
grid on
% Plot F3
subplot(2,2,4)
plot(t2*180/pi, F3)
xlabel('Theta2, deg')
ylabel('F3, N-m')
axis tight
grid on
