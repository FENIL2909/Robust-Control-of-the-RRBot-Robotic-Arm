clc; clear; close;
syms theta1 theta2 theta1_dot theta2_dot theta1_ddot theta2_ddot tau1 tau2 g 'real'
syms t 'real'
visualization = true;

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
%% Virtual Control input design

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
B = [0 0; 0 0; 1 0; 0 1];

lambda = [-3 -3 -4 -4];
fprintf("**************************************************************************************************\n")
fprintf("********** K Gains **********")

K = place(A,B,lambda)
Kp = K(:,1:2)             %[12 0; 0 12]
Kd = K(:,3:4)
O = [0 0; 0 0];
Acl = [O eye(2); -Kp, -Kd];
Q = eye(4)*20;
P = lyap(Acl',Q)
rho = 3.25;
phi = 0.075;

T = 10;
theta1_initial = 200;
theta2_initial = 125;

y0 = [deg2rad(theta1_initial), deg2rad(theta2_initial) 0, 0];

[time,y] = ode45(@ode_rrbot,[0,T],y0);

theta1_desired = (pi*time.^3)/500 - (3*pi*time.^2)/100 - time/18014398509481984 + pi;
theta2_desired = (pi*time.^3)/1000 - (3*pi*time.^2)/200 - time/36028797018963968 + pi/2;

theta1_dot_desired = (3*pi*time.^2)/500 - (3*pi*time)/50 - 1/18014398509481984;
theta2_dot_desired = (3*pi*time.^2)/1000 - (3*pi*time)/100 - 1/36028797018963968;

theta1_ddot_desired = (3*pi*time)/250 - (3*pi)/50;
theta2_ddot_desired = (3*pi*time)/500 - (3*pi)/100;

a_hat = I1_hat + I2_hat + M1_hat*r1^2 + M2_hat*(L1^2 + r2^2);
b_hat = M2_hat*L1*r2;
d_hat = I2_hat + M2_hat*r2^2;


% Reconstructing the Control Inputs

for i = 1:size(time)
    
    Mmat_hat= [a_hat+2*b_hat*cos(y(i,2)), d_hat+b_hat*cos(y(i,2)); d_hat+b_hat*cos(y(i,2)), d_hat];
    Cmat_hat= [-b_hat*sin(y(i,2))*y(i,4), -b_hat*sin(y(i,2))*(y(i,1) + y(i,2)); b_hat*sin(y(i,2))*y(i,3),0];
    Gmat_hat= [-M1_hat*g*r1*sin(y(i,1))-M2_hat*g*(L1*sin(y(i,1))+r2*sin(y(i,1)+y(i,2))); -M2_hat*g*r2*sin(y(i,1)+y(i,2))];

    theta1_desired_1 = (pi*time(i)^3)/500 - (3*pi*time(i)^2)/100 - time(i)/18014398509481984 + pi;
    theta2_desired_1 = (pi*time(i)^3)/1000 - (3*pi*time(i)^2)/200 - time(i)/36028797018963968 + pi/2;
    theta1_dot_desired_1 = (3*pi*time(i)^2)/500 - (3*pi*time(i))/50 - 1/18014398509481984;
    theta2_dot_desired_1 = (3*pi*time(i)^2)/1000 - (3*pi*time(i))/100 - 1/36028797018963968;
    theta1_ddot_desired_1 = (3*pi*time(i))/250 - (3*pi)/50;
    theta2_ddot_desired_1 = (3*pi*time(i))/500 - (3*pi)/100;

    feed_foward_input = [theta1_ddot_desired_1; theta2_ddot_desired_1];
    e = [y(i,1) - theta1_desired_1; y(i,2) - theta2_desired_1];
    e_dot = [y(i,3) - theta1_dot_desired_1; y(i,4) - theta2_dot_desired_1];

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

%     Vr = [0;0];

    V = feed_foward_input - Kp*e - Kd*e_dot + Vr;
    U = Mmat_hat*V + Cmat_hat*[y(i,3); y(i,4)] + Gmat_hat;
    U1(i) = U(1);
    U2(i) = U(2);
end


%% plots
figure
hold on
subplot(2,2,1)
plot(time,y(:,1))
hold on
plot(time, theta1_desired,'Color','red','LineStyle','--')
xlabel('Time step')
ylabel('rad')
title('theta1')

subplot(2,2,2)
plot(time,y(:,2))
hold on
plot(time, theta2_desired,'Color','red','LineStyle','--')
xlabel('Time step')
ylabel('rad')
title('theta2')

subplot(2,2,3)
plot(time,y(:,3))
hold on
plot(time, theta1_dot_desired,'Color','red','LineStyle','--')
xlabel('Time step')
ylabel('rad/s')
title('theta1-dot')

subplot(2,2,4)
plot(time,y(:,4))
hold on
plot(time, theta2_dot_desired,'Color','red','LineStyle','--')
xlabel('Time step')
ylabel('rad/s')
title('theta2-dot')
hold off

figure
hold on
subplot(2,1,1)
plot(time,U1)
xlabel('Time step')
ylabel('Nm')
title('tau1')

subplot(2,1,2)
plot(time,U2)
xlabel('Time step')
ylabel('Nm')
title('tau2')
hold off
pause(10)
close all

%% visualization
figure
if(visualization)

    theta1_desired_plot = (pi*time.^3)/500 - (3*pi*time.^2)/100 - time/18014398509481984 + pi;
    theta2_desired_plot = (pi*time.^3)/1000 - (3*pi*time.^2)/200 - time/36028797018963968 + pi/2;
   
    l1 = 1;
    l2 = 1;

    x1_desired_plot = l1*sin(theta1_desired_plot);
    y1_desired_plot = l1*cos(theta1_desired_plot);

    x2_desired_plot = l1*sin(theta1_desired_plot) + l2*sin(theta1_desired_plot+theta2_desired_plot);
    y2_desired_plot = l1*cos(theta1_desired_plot) + l2*cos(theta1_desired_plot+theta2_desired_plot);
    

    x1_pos= l1*sin(y(:,1));
    x2_pos= l1*sin(y(:,1)) + l2*sin(y(:,1)+y(:,2));
    y1_pos= l1*cos(y(:,1));
    y2_pos= l1*cos(y(:,1)) + l2*cos(y(:,1)+y(:,2));

    for i=1:size(y)
        plot(x1_desired_plot, y1_desired_plot,'Color','red','LineStyle','--')
        hold on
        plot(x2_desired_plot, y2_desired_plot,'Color','red','LineStyle','--')
        hold on
        plot([0 x1_pos(i) x2_pos(i)],[0 y1_pos(i) y2_pos(i)],'blue', 'LineWidth',4.0)
        hold on
        plot(x2_pos(1:i),y2_pos(1:i),'blue')
        hold on
        plot(x1_pos(1:i),y1_pos(1:i),'blue')
        hold on
        plot(x1_pos(i),y1_pos(i),'.','MarkerSize',24.0)
        xlim([-2.25 2.25])
        ylim([-2.25 2.25])
        pause(0.000075)
        hold off
    end
end
