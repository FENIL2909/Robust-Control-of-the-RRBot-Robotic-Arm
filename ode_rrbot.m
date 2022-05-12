%% ode_rrbot
function dX = ode_rrbot(t,X)

boundary_layer = true;
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

dX = zeros(4,1);
X = num2cell(X);
[theta1, theta2, theta1_dot, theta2_dot] = deal(X{:});

a_hat = I1_hat + I2_hat + M1_hat*r1^2 + M2_hat*(L1^2 + r2^2);
b_hat = M2_hat*L1*r2;
d_hat = I2_hat + M2_hat*r2^2;

Mmat_hat= [a_hat+2*b_hat*cos(theta2), d_hat+b_hat*cos(theta2); d_hat+b_hat*cos(theta2), d_hat];
Cmat_hat= [-b_hat*sin(theta2)*theta2_dot, -b_hat*sin(theta2)*(theta1_dot + theta2_dot); b_hat*sin(theta2)*theta1_dot,0];
Gmat_hat= [-M1_hat*g*r1*sin(theta1)-M2_hat*g*(L1*sin(theta1)+r2*sin(theta1+theta2)); -M2_hat*g*r2*sin(theta1+theta2)];


theta1_desired = (pi*t^3)/500 - (3*pi*t^2)/100 - t/18014398509481984 + pi;
theta2_desired = (pi*t^3)/1000 - (3*pi*t^2)/200 - t/36028797018963968 + pi/2;

theta1_dot_desired = (3*pi*t^2)/500 - (3*pi*t)/50 - 1/18014398509481984;
theta2_dot_desired = (3*pi*t^2)/1000 - (3*pi*t)/100 - 1/36028797018963968;

theta1_ddot_desired = (3*pi*t)/250 - (3*pi)/50;
theta2_ddot_desired = (3*pi*t)/500 - (3*pi)/100;

feed_foward_input = [theta1_ddot_desired; theta2_ddot_desired];
e = [theta1 - theta1_desired; theta2 - theta2_desired];
e_dot = [theta1_dot - theta1_dot_desired; theta2_dot - theta2_dot_desired];
%% Virtual Control input design

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
B = [0 0; 0 0; 1 0; 0 1];

lambda = [-3 -3 -4 -4];

K = place(A,B,lambda);
Kp = K(:,1:2);              %[12 0; 0 12]
Kd = K(:,3:4);
O = [0 0; 0 0];
Acl = [O eye(2); -Kp, -Kd] ;
Q = eye(4)*20;
P = lyap(Acl',Q);
rho = 3.25;
phi = 0.075;

G = [e; e_dot];

if boundary_layer
    if norm(B'*P*G) > phi
        Vr = -(rho*(B'*P*G))/norm(B'*P*G);
    else
        Vr = -(rho*(B'*P*G))/phi;
    end
else
    if norm(B'*P*G) ~= 0
        Vr = -(rho*(B'*P*G))/norm(B'*P*G);
    else
        Vr = [0;0];
    end
end

% Vr = [0;0];

V = feed_foward_input - Kp*e - Kd*e_dot + Vr;

U = Mmat_hat*V + Cmat_hat*[theta1_dot; theta2_dot] + Gmat_hat;

tau1 = U(1);
tau2 = U(2);

dX(1) = theta1_dot;
dX(2) = theta2_dot;
dX(3) = (I2*tau1 - I2*tau2 + M2*r2^2*tau1 - M2*r2^2*tau2 + L1*M2^2*g*r2^2*sin(theta1) + I2*L1*M2*g*sin(theta1) + I2*M1*g*r1*sin(theta1) - L1*M2*r2*tau2*cos(theta2) + L1*M2^2*r2^3*theta1_dot^2*sin(theta2) + L1*M2^2*r2^3*theta2_dot^2*sin(theta2) + L1^2*M2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - L1*M2^2*g*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*L1*M2*r2*theta1_dot^2*sin(theta2) + I2*L1*M2*r2*theta2_dot^2*sin(theta2) + M1*M2*g*r1*r2^2*sin(theta1) + 2*L1*M2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*I2*L1*M2*r2*theta1_dot*theta2_dot*sin(theta2))/(- L1^2*M2^2*r2^2*cos(theta2)^2 + L1^2*M2^2*r2^2 + I2*L1^2*M2 + M1*M2*r1^2*r2^2 + I1*M2*r2^2 + I2*M1*r1^2 + I1*I2);
dX(4) = -(I2*tau1 - I1*tau2 - I2*tau2 - L1^2*M2*tau2 - M1*r1^2*tau2 + M2*r2^2*tau1 - M2*r2^2*tau2 - L1^2*M2^2*g*r2*sin(theta1 + theta2) + L1*M2^2*g*r2^2*sin(theta1) - I1*M2*g*r2*sin(theta1 + theta2) + I2*L1*M2*g*sin(theta1) + I2*M1*g*r1*sin(theta1) + L1*M2*r2*tau1*cos(theta2) - 2*L1*M2*r2*tau2*cos(theta2) + L1*M2^2*r2^3*theta1_dot^2*sin(theta2) + L1^3*M2^2*r2*theta1_dot^2*sin(theta2) + L1*M2^2*r2^3*theta2_dot^2*sin(theta2) + 2*L1^2*M2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + L1^2*M2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - L1*M2^2*g*r2^2*sin(theta1 + theta2)*cos(theta2) + L1^2*M2^2*g*r2*cos(theta2)*sin(theta1) - M1*M2*g*r1^2*r2*sin(theta1 + theta2) + I1*L1*M2*r2*theta1_dot^2*sin(theta2) + I2*L1*M2*r2*theta1_dot^2*sin(theta2) + I2*L1*M2*r2*theta2_dot^2*sin(theta2) + M1*M2*g*r1*r2^2*sin(theta1) + 2*L1*M2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*L1^2*M2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + L1*M1*M2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*L1*M2*r2*theta1_dot*theta2_dot*sin(theta2) + L1*M1*M2*g*r1*r2*cos(theta2)*sin(theta1))/(- L1^2*M2^2*r2^2*cos(theta2)^2 + L1^2*M2^2*r2^2 + I2*L1^2*M2 + M1*M2*r1^2*r2^2 + I1*M2*r2^2 + I2*M1*r1^2 + I1*I2);
end