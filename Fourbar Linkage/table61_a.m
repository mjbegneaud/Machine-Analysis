% Author: Sally McInerny, Ph.D. 
% Modified: Matthew J. Begneaud
%This routine solves for the unknown thetas, theta3 and theta4
% for a GCRR for a full range rotation of the crank
clear
clc

global a b c d t2 w2 %This allows me to set the values of these
                                % variables in functions that I call
                                % without having to pass them as input
                                % parameters
t2start=0; %This is the value of theta2 in degrees for which I have 
t2stop=360; 
% a graphical solution that provides estimates of the
                      % corresponding theta3 and theta4
t3est=106.6;  %Measured theta3 in degrees
t4est=131.8; %Measured theta4 in degrees
tunknownguess=[t3est; t4est]*pi/180; % vector of the two initial guesses of 
                    % theta3 and theta4,in radians
a=2;
b=7;
c=9;
d=6;

w2=10;
a2=0



% The first calculation uses the specified starting value of theta2 and
% provides initial estimates of theta3 and theta4 from the graphical
% solution
%
% After that, I am going to increment theta2 by five degrees, and use the
% values of theta3 and theta4 from the previous calculation as the
% initial estimates for the new ones
% 
%  I'm going to repeat (iterate) this procedure until I've returned back to
t2=t2start*pi/180;
n=floor(((t2stop+5)-t2start)/5); %in case this is non-integer, the value is truncated at 
                               % whole number to make it an integer
i=1;
unknownthetas=zeros(n,2); %preallocates memory space for the matrix
Fval=zeros(n,2); 

% w2val = ones(n,1);
% w2val = w2*w2val;

while (1) % Will continue until the condition to break 
                  %  out of the while loop is met
[unknownthetas(i,:),Fval(i,:)]=fsolve(@thetwoeqtns,tunknownguess);
t2=t2+5*pi/180;
tunknownguess=unknownthetas(i,:);
i=i+1;
if i>n
    break
end
end

t2=[t2start:5:t2stop]'*pi/180; %Generates the full vector of 
                        %  theta2 values.  This could've been done along
                        %  the way, but it would've required that the
                        %  current value of t2 be passed as an additional
                        %  parameter, instead of relying on making it a
                        %  global variable
t3 = unknownthetas(:,1);
t4 = unknownthetas(:,2);

omegas = zeros(2,n);
for i = 1:n
    A = [-b*sin(t3(i)), c*sin(t4(i));...
        b*cos(t3(i)), -c*cos(t4(i))];
    rhs = a*w2*[sin(t2(i)), -cos(t2(i))]';
    omegas(:,i) = A\rhs;
end

w3 = omegas(1,:)';
w4 = omegas(2,:)';

alphas = zeros(2,n);
for i = 1:n
    B = [-b*sin(t3(i)), c*sin(t4(i));...
        b*cos(t3(i)), -c*cos(t4(i))];
    arhs = a*a2*[sin(t2(i)), -cos(t2(i))]' + a*(w2^2)*[cos(t2(i)), sin(t2(i))]'...
        + b*(w3(i)^2)*[cos(t3(i)), sin(t3(i))]' - c*(w4(i)^2)*[cos(t4(i)), sin(t4(i))]';
    alphas(:, i) = B\arhs;
end

a3 = alphas(1,:)';
a4 = alphas(2,:)';

%Plot w3 vs theta2 and w4 vs theta2 
subplot(2, 2, 1)
plot(t2*180/pi, w3)
xlabel('Theta2, deg')
ylabel('Omega3, rad/s')
axis tight
grid on
subplot(2,2,2)
plot(t2*180/pi, w4)
xlabel('Theta2, deg')
ylabel('Omega4, rad/s')
axis tight
grid on 

%Plot a3 vs theta2 and a4 vs theta2
subplot(2,2,3)
plot(t2*180/pi, a3)
xlabel('Theta2, deg')
ylabel('Alpha3, rad/s^2')
axis tight
grid on
subplot(2,2,4)
plot(t2*180/pi, a4)
xlabel('Theta2, deg')
ylabel('Alpha4, rad/s^2')
axis tight
grid on




%Plot theta3 versus theta2 and theta4 versus theta2 in separate plots
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

