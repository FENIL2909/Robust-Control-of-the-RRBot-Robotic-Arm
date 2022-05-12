%% System Parameters

M1 = 1; %Kg
M2 = 1; %Kg
L1 = 1; %m
L2 = 1; %m
r1 = 0.45; %m
r2 = 0.45; %m
I1 = 0.084; %Kg.m2
I2 = 0.084; %Kg.m2
g = 9.81; %m/s2

M1_hat = 0.75; %Kg
M2_hat = 0.75; %Kg
I1_hat = 0.063; %Kg.m2
I2_hat = 0.063;%Kg.m2

syms theta1 theta2 theta1_dot theta2_dot theta1_ddot theta2_ddot tau1 tau2 g 'real'
syms M1 M2 M1_hat M2_hat I1 I2 I1_hat I2_hat L1 L2 r1 r2 'real'
syms x1_dot y1_dot x2_dot y2_dot 'real'
syms t 'real'

%% State Space Representation
X = sym('X', [4,1]);
X(1) = theta1;
X(2) = theta2;
X(3) = theta1_dot;
X(4) = theta2_dot;

%% Dynamic Equations

x1_dot = (theta1_dot)*(r1)*(cos(theta1));
y1_dot = -(theta1_dot)*(r1)*(sin(theta1));

x2_dot = (theta1_dot)*(L1)*(cos(theta1)) + (theta1_dot + theta2_dot)*(r2)*(cos(theta1 + theta2));
y2_dot = -(theta1_dot)*(L1)*(sin(theta1)) - (theta1_dot + theta2_dot)*(r2)*(sin(theta1 + theta2));

K1 = (1/2)*(I1)*(theta1_dot*theta1_dot) + (1/2)*(M1)*((x1_dot*x1_dot) + (y1_dot*y1_dot));
K2 = (1/2)*(I2)*((theta2_dot + theta1_dot)*(theta2_dot + theta1_dot)) + (1/2)*(M2)*((x2_dot*x2_dot) + (y2_dot*y2_dot));

P1 = M1*g*r1*cos(theta1);
P2 = M2*g*(L1*(cos(theta1)) + r2*(cos(theta1 + theta2)));
L = K1 + K2 - P1 - P2;

u = [tau1;tau2];
q = [theta1;theta2];
dq = [theta1_dot; theta2_dot];
ddq = [theta1_ddot; theta2_ddot];

DL_Dq = gradient(L,q);  % used gradient instead of jacobian to keep matrix size consistent
DL_Ddq = gradient(L,dq); % used gradient instead of jacobian to keep matrix size consistent
dDL_dtDdq = jacobian(DL_Ddq,[q;dq])*[dq;ddq];

EOM = simplify(dDL_dtDdq - DL_Dq -u);
%disp(EOM)

EOM_numerical = subs(EOM,[M1,M2,L1,L2,r1,r2,I1,I2,g],[1,1,1,1,0.45,0.45,0.084,0.084,9.81]);


%% Feedback Linearization

a = I1 + I2 + M1*r1^2 + M2*(L1^2 + r2^2);
b = M2*L1*r2;
d = I2 + M2*r2^2;

Mmat= [a+2*b*cos(theta2), d+b*cos(theta2); d+b*cos(theta2), d];
Cmat= [-b*sin(theta2)*theta2_dot, -b*sin(theta2)*(theta1_dot + theta2_dot); b*sin(theta2)*theta1_dot,0];
Gmat= [-M1*g*r1*sin(theta1)-M2*g*(L1*sin(theta1)+r2*sin(theta1+theta2)); -M2*g*r2*sin(theta1+theta2)];

a_hat = I1_hat + I2_hat + M1_hat*r1^2 + M2_hat*(L1^2 + r2^2);
b_hat = M2_hat*L1*r2;
d_hat = I2_hat + M2_hat*r2^2;

Mmat_hat= [a_hat+2*b_hat*cos(theta2), d_hat+b_hat*cos(theta2); d_hat+b_hat*cos(theta2), d_hat];
Cmat_hat= [-b_hat*sin(theta2)*theta2_dot, -b_hat*sin(theta2)*(theta1_dot + theta2_dot); b_hat*sin(theta2)*theta1_dot,0];
Gmat_hat= [-M1_hat*g*r1*sin(theta1)-M2_hat*g*(L1*sin(theta1)+r2*sin(theta1+theta2)); -M2_hat*g*r2*sin(theta1+theta2)];

fprintf("**************************************************************************************************\n")
fprintf("********** EOM Actual **********")
EOM_FL = Mmat*[theta1_ddot; theta2_ddot] + Cmat*[theta1_dot; theta2_dot] + Gmat

fprintf("**************************************************************************************************\n")
fprintf("**********EOM with Nominal Parameters **********")
EOM_FL_hat = Mmat_hat*[theta1_ddot; theta2_ddot] + Cmat_hat*[theta1_dot; theta2_dot] + Gmat_hat




