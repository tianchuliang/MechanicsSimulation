
% Define system of equations
% Turns out that if m1 = m2, then the equations of motion are independent
% of mass! These equations only work for m1 = m2 = m and L1 = L2 = L.
function [ df ] = doublepend( t,q,L,g)
df = zeros(4,1);
df(1) = q(3); % dq1 = q1dot
df(2) = q(4); % dq2 = q2dot
df(3) =(3*(9*g*sin(q(1)) + 3*g*sin(q(1)-2*q(2)) + ...
    2*L*(2*q(4)^2 + 3*q(3)^2*cos(q(1) - q(2)))*sin(q(1) -...
    q(2))))/(L*(-23 + 9*cos(2*(q(1) - q(2))))); %dq1dot
df(4) = -((3*(16*L*q(3)^2*sin(q(1) - q(2)) + ...
    3*L*q(4)^2*sin(2*(q(1) - q(2))) + 9*g*sin(2*q(1) -...
    q(2)) - 7*g*sin(q(2))))/(L*(-23 + 9*cos(2*(q(1) - q(2)))))); % dq2dot
end

