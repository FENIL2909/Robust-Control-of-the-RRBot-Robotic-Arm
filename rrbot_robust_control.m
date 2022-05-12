clear; close; clc;
% ROS Setup
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
rosinit;
j1_effort = rospublisher('/rrbot/joint1_effort_controller/command');
j2_effort = rospublisher('/rrbot/joint2_effort_controller/command');
JointStates = rossubscriber('/rrbot/joint_states');
tau1 = rosmessage(j1_effort);
tau2 = rosmessage(j2_effort);
tau1.Data = 0;
tau2.Data = 0;
send(j1_effort,tau1);
send(j2_effort,tau2);
client = rossvcclient('/gazebo/set_model_configuration');
req = rosmessage(client);
req.ModelName = 'rrbot';
req.UrdfParamName = 'robot_description';
req.JointNames = {'joint1','joint2'};
req.JointPositions = [deg2rad(200), deg2rad(125)];
resp = call(client,req,'Timeout',3);

i = 1;

tic;
t = 0;
while(t < 10)
t = toc;
% read the joint states
jointData = receive(JointStates);
% inspect the "jointData" variable in MATLAB to get familiar with its structure
% design your state feedback controller in the following

theta1 = wrapTo2Pi(jointData.Position(1));
theta2 = wrapTo2Pi(jointData.Position(2));
theta1_dot = jointData.Velocity(1);
theta2_dot = jointData.Velocity(2);

X = [theta1;theta2;theta1_dot;theta2_dot];

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

V = feed_foward_input - Kp*e - Kd*e_dot + Vr;

U = Mmat_hat*V + Cmat_hat*[theta1_dot; theta2_dot] + Gmat_hat;

tau1.Data = U(1);
tau2.Data = U(2);
send(j1_effort,tau1);
send(j2_effort,tau2);

% you can sample data here to be plotted at the end
X1(i) = theta1;
X2(i) = theta2;
X3(i) = theta1_dot;
X4(i) = theta2_dot;
TAU1(i) = U(1);
TAU2(i) = U(2);

time(i) = t;
i = i+1;

end
tau1.Data = 0;
tau2.Data = 0;
send(j1_effort,tau1);
send(j2_effort,tau2);
% disconnect from roscore
rosshutdown;

Theta1_desired = (pi*time.^3)/500 - (3*pi*time.^2)/100 - time/18014398509481984 + pi;
Theta2_desired = (pi*time.^3)/1000 - (3*pi*time.^2)/200 - time/36028797018963968 + pi/2;

Theta1_dot_desired = (3*pi*time.^2)/500 - (3*pi*time)/50 - 1/18014398509481984;
Theta2_dot_desired = (3*pi*time.^2)/1000 - (3*pi*time)/100 - 1/36028797018963968;


%% plots
figure
hold on
subplot(2,2,1)
plot(time,X1)
hold on
plot(time, Theta1_desired,'Color','red','LineStyle','--')
xlabel('Time step')
ylabel('rad')
title('theta1')

subplot(2,2,2)
plot(time,X2)
hold on
plot(time, Theta2_desired,'Color','red','LineStyle','--')
xlabel('Time step')
ylabel('rad')
title('theta2')

subplot(2,2,3)
plot(time,X3)
hold on
plot(time, Theta1_dot_desired,'Color','red','LineStyle','--')
xlabel('Time step')
ylabel('rad/s')
title('theta1-dot')

subplot(2,2,4)
plot(time,X4)
hold on
plot(time, Theta2_dot_desired,'Color','red','LineStyle','--')
xlabel('Time step')
ylabel('rad/s')
title('theta2-dot')
hold off

figure
hold on
subplot(2,1,1)
plot(time,TAU1)
xlabel('Time step')
ylabel('Nm')
title('tau1')

subplot(2,1,2)
plot(time,TAU2)
xlabel('Time step')
ylabel('Nm')
title('tau2')
hold off