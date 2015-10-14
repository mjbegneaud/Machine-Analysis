% Author: Sally McInerny, Ph.D. 
% Modified: Matthew J. Begneaud
%This routine solves for the unknown thetas, theta3, theta4, omega3,
%omega4, alpha3, and alpha 4 for a fourbar linkage
clear
clc

% This allows user to set the values of these variables in functions that 
% are called without having to pass them as input parameters
global a b c d t2 w2 a2

% Angular range of motion (in deg) for link2 (Theta2 range)
t2start=-70; 
t2stop=70; 

% Initial guess (in deg) for Theta3 and Theta4 (vector of guesses included)
t3est=47.202;
t4est=191.884;
tunknownguess=[t3est; t4est]*pi/180; 

% Link lengths
a=10;
b=10;
c=10;
d=20;

% Omega2 and Alpha2
% w2=10;
% a2=0;



% The first calculation uses the specified starting value of theta2 and
% provides initial estimates of theta3 and theta4 from the graphical
% solution
%
% After that, I am going to increment theta2 by five degrees, and use the
% values of theta3 and theta4 from the previous calculation as the
% initial estimates for the new ones

% Setup for step size of Theta2
t2=t2start*pi/180;
n=floor(((t2stop+5)-t2start)/5); 

% Start counter
i=1;

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
% omegas = zeros(2,n);
% for i = 1:n
%     A = [-b*sin(t3(i)), c*sin(t4(i));...
%         b*cos(t3(i)), -c*cos(t4(i))];
%     rhs = a*w2*[sin(t2(i)), -cos(t2(i))]';
%     omegas(:,i) = A\rhs;
% end
% % Generate full vector of Omega3 and Omega4 values
% w3 = omegas(1,:)';
% w4 = omegas(2,:)';
% 
% 
% % Main loop for Alphas
% alphas = zeros(2,n);
% for i = 1:n
%     B = [-b*sin(t3(i)), c*sin(t4(i));...
%         b*cos(t3(i)), -c*cos(t4(i))];
%     arhs = a*a2*[sin(t2(i)), -cos(t2(i))]' + a*(w2^2)*[cos(t2(i)), sin(t2(i))]'...
%         + b*(w3(i)^2)*[cos(t3(i)), sin(t3(i))]' - c*(w4(i)^2)*[cos(t4(i)), sin(t4(i))]';
%     alphas(:, i) = B\arhs;
% end
% % Generate full vector of Alpha3 and Alpha4 values
% a3 = alphas(1,:)';
% a4 = alphas(2,:)';


% Plot theta3 versus theta2 and theta4 versus theta2 in separate plots
subplot(2,1,1)
plot(t2*180/pi,unknownthetas(:,1)*180/pi)
axis tight
grid on
xlabel('Theta2, deg')
ylabel('Theta3, deg')
subplot(2,1,2)
plot(t2*180/pi,unknownthetas(:,2)*180/pi)
axis tight
grid on
xlabel('Theta2, deg')
ylabel('Theta4,deg')


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




