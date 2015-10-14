% Matthew J. Begneaud, 9/30/15
% Script to analyze the output of the rocker-link in a crank-slider-rocker 

clear
clc
    
global a d t2 w2 

a = 5;
d = 20;
w2 = 10;

% Theta2 range in degrees
t2start = 0;
t2stop = 360;
n=floor(((t2stop+5)-t2start)/5); 

t2 = (t2start:5:t2stop)*pi/180;




% %%%%%%%%%%%%%%%%%%%%%%%%% Position %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t4 = atan2(d+a*sin(t2), a*cos(t2));

c = zeros(1,n);
for i = 1:n
    c (:,i) = (d+a*sin(t2(:,i)))/(sin(t4(:,i)));
end




% %%%%%%%%%%%%%%%%%%%%%%%%%% Velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
velocity = zeros(2,n);
for i = 1:n
    A = [cos(t4(i)) , -c(i).*sin(t4(i));...
        sin(t4(i)) , c(i).*cos(t4(i))];
    vrhs = a*w2*[-sin(t2(i)) ; cos(t2(i))];
    
    velocity(:,i) =  A\vrhs;
end
cdot = velocity(1,:);
w4 = velocity(2,:);




% %%%%%%%%%%%%%%%%%%%%%%%%%%% Acceleration %%%%%%%%%%%%%%%%%%%%%%%%%%
acceleration = zeros(2,n);
for i = 1:n
    B = [cos(t4(i)) , -c(i)*sin(t4(i));...
        sin(t4(i)) , c(i)*cos(t4(i))];
    arhs = -a*w2^2*[cos(t2(i)); sin(t2(i))] + c(i)*w4(i)^2*[cos(t4(i)); sin(t4(i))] + 2*cdot(i)*w4(i)*[sin(t4(i)); -cos(t4(i))];
    
    acceleration (:,i) = B\arhs;
end
cddot = acceleration(1,:);
a4 = acceleration(2,:);




% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot c and t4 vs t2
subplot(3,2,1)
plot(t2*180/pi, t4*180/pi)
xlabel('Theta2, deg')
ylabel('Theta4, deg')
axis tight
grid on
subplot(3,2,2)
plot(t2*180/pi, c)
xlabel('Theta2, deg')
ylabel('C, m')
axis tight
grid on 

% Plot cdot and w4 vs t2
subplot(3,2,3)
plot(t2*180/pi, w4*180/pi)
xlabel('Theta2, deg')
ylabel('Omega4, deg/s')
axis tight
grid on
subplot(3,2,4)
plot(t2*180/pi, cdot)
xlabel('Theta2, deg')
ylabel('Cdot, m/s')
axis tight
grid on

% Plot cddot and a4 vs t2
subplot(3,2,5)
plot(t2*180/pi, a4*180/pi)
xlabel('Theta2, deg')
ylabel('Alpha4, deg/s^2')
axis tight
grid on
subplot(3,2,6)
plot(t2*180/pi, cddot)
xlabel('Theta2, deg')
ylabel('Cddot, m/s^2')
axis tight
grid on