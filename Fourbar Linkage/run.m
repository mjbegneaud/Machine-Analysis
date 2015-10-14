clear
clc
%Example 4-1 pg 190 of Norton, 5th ed
t3est=15*pi/180;
t4est=60*pi/180;
theta3_4_guess=[t3est; t4est];
global a b c d t2
t2=40*pi/180;
a=40;
b=120;
c=80;
d=100;

[theta3_4,Fval]=fsolve('thetwoeqtns',theta3_4_guess)
% Print calculated values to command screen in degrees
%  represent the answer to two decimal places
round(100*theta3_4'*180/pi)/100
% Note that I've taken the transpose of the column vector
%  so that the two values are appear side by side in a row 
