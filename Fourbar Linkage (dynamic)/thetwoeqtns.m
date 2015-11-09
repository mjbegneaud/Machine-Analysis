function F=thetwoeqtns(theta3_4)
global a b c d t2
% The fsolve routine expects the n-equations to be
%  in nx1 column vector form
F=[d+c*cos(theta3_4(2))-a*cos(t2)-b*cos(theta3_4(1)); ...
    c*sin(theta3_4(2))-a*sin(t2)-b*sin(theta3_4(1))];
