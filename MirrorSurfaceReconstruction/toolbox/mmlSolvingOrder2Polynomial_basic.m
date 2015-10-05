function [sol_final coefmat]= mmlSolvingOrder2Polynomial_basic(x,y,w1x,w2x,w3x,w1y,w2y,w3y,w1,w2,w3)

% define the symbols
% w = (w1;w2;w3);
% w1x: partial derivative of w1 over x

sol_final = [];


% syms w1x w2x w3x w1 w2 w3 w1y w2y w3y x y s

% corresponding to #1
poly1 = (x*x+y*y+1);

% corresponding to #3
poly3 = 2*w3+2*x*w1+2*y*w2;

% corresponding to A
tmpA = w2x*poly1^2 - w1y*poly1^2+w2x*x^2*poly1^2-w1y*y^2*poly1^2+w2*x*poly1-w1*y*poly1-w2*x*poly1^2+w2*x^3*poly1+w3y*x*poly1^2+w1*y*poly1^2-w1*y^3*poly1-w3x*y*poly1^2-w1*x^2*y*poly1+w2*x*y^2*poly1-w1x*x*y*poly1^2+w2y*x*y*poly1^2;

% corresponding to a
tmpa = w2x*poly1^(3/2)-w1y*poly1^(3/2)+w2x*x^2*poly1^(3/2)-w1y*y^2*poly1^(3/2)+w3y*x*poly1^(3/2)-w3x*y*poly1^(3/2)-w1x*x*y*poly1^(3/2)+w2y*x*y*poly1^(3/2);

% corresponding to B
tmpB = w1*y*poly3 - w2*x*poly3+w2*w3x*poly1-w1*w3y*poly1-w2*x^3*poly3+w1*y^3*poly3-w2x*poly1*poly3+w1y*poly1*poly3+w2*x*poly1*poly3-w3y*x*poly1*poly3-w1*y*poly1*poly3+w3x*y*poly1*poly3+w2*w1x*x^3*poly1+w2*w3x*x^2*poly1...
    +w3*w1y*x^2*poly1-w3*w2x*y^2*poly1-w1*w2y*y^3*poly1-w1*w3y*y^2*poly1-w2x*x^2*poly1*poly3+w1y*y^2*poly1*poly3+w2*w1x*x*poly1-w1*w1y*x*poly1+w3*w3y*x*poly1+w2*w2x*y*poly1-w3*w3x*y*poly1-w1*w2y*y*poly1+w1*x^2*y*poly3...
    -w2*x*y^2*poly3-w1*w3x*x*y*poly1-w3*w1x*x*y*poly1+w2*w3y*x*y*poly1+w3*w2y*x*y*poly1+w1x*x*y*poly1*poly3-w2y*x*y*poly1*poly3-w1*w1x*x^2*y*poly1-w1*w2x*x*y^2*poly1+w2*w2x*x^2*y*poly1-w1*w1y*x*y^2*poly1+w2*w1y*x^2*y*poly1...
    +w2*w2y*x*y^2*poly1;

% corresponding to b
tmpb = w2*w3x*poly1^(3/2)-w3*w2x*poly1^(3/2)-w1*w3y*poly1^(3/2)+w3*w1y*poly1^(3/2)-w1*w2x*x*poly1^(3/2)+w2*w1x*x*poly1^(3/2)-w1*w2y*y*poly1^(3/2)+w2*w1y*y*poly1^(3/2);

% the coefficient of s in #2
tmpQ = -poly3;

% corresponding to #4
poly4 = w1^2 + w2^2 + w3^2;

% corresponding to the constant in #2
tmpT = poly4;

% corresponding to the coefficient of s^2 in #2
tmpP = poly1;

% corresponding to C
tmpC = (w2*x -w1*y)*poly4+(w2*x^3 - w1*y^3)*poly4+(-w2^2*w2x+w1^2*w1y)*poly1+(w2x-w1y)*poly1*poly4+(-w1*x+w2*y)*x*y*poly4+(-w3y*x+w3x*y)*w3^2*poly1+(-w2*x+w3y*x+w1*y-w3x*y)*poly1*poly4+(-w2^2*w2x*x^2+w1^2*w1y*y^2)*poly1+...
    (w2x*x^2-w1y*y^2)*poly1*poly4+(-w1*w1x-w3*w3x)*w2*poly1+(w2*w2y+w3*w3y)*w1*poly1+(-w1*w1y-w2*w2y)*w3*x*poly1+(w1*w1x+w2*w2x)*w3*y*poly1+(-w1*w1x-w3*w3x)*w2*x^2*poly1+(w2*w2y+w3*w3y)*w1*y^2*poly1+(w1^2*w1x-w2^2*w2y)*x*y*poly1+...
    (-w1x+w2y)*x*y*poly1*poly4+(w2*w2x+w3*w3x)*w1*x*y*poly1+(-w1*w1y-w3*w3y)*w2*x*y*poly1;

% coefficient for s^2
testorder2 = tmpB^2+2*tmpA*tmpC-tmpb^2*tmpP-2*tmpa*tmpb*tmpQ-tmpa^2*tmpT;

% coefficient for s
testorder1 = 2*tmpB*tmpC - tmpb^2*tmpQ - 2*tmpa*tmpb*tmpT;

% constant
testorder0 = tmpC^2 - tmpb^2*tmpT;

coefmat = [testorder2 testorder1 testorder0];

% solutions: solution one and solution two
sol_s_1 = (-testorder1+sqrt(testorder1^2 - 4*testorder2*testorder0))/(2*testorder2);
sol_s_2 = (-testorder1-sqrt(testorder1^2 - 4*testorder2*testorder0))/(2*testorder2);

tmpval_s1 = -(tmpA*sol_s_1^2+tmpB*sol_s_1+tmpC)/(tmpa*sol_s_1+tmpb);
tmpval_s2 = -(tmpA*sol_s_2^2+tmpB*sol_s_2+tmpC)/(tmpa*sol_s_2+tmpb);


% basic test
if sol_s_1 >0 && tmpval_s1 >0 && isreal(sol_s_1)==1,
    sol_final = [sol_final sol_s_1];
 
end

if sol_s_2 > 0 && tmpval_s2 > 0 && isreal(sol_s_2)==1,
    sol_final = [sol_final sol_s_2];
end


end



