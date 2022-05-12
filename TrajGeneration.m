%% Trajectory Generation
syms a01 a11 a21 a31 a02 a12 a22 a32 t 'real'

% Traj_Coeff1 = [a01 a11 a21 a31]';
% Traj_Coeff2 = [a02 a12 a22 a32]';

theta1_initial = deg2rad(180);
theta2_initial = deg2rad(90);
theta1_final = 0;
theta2_final = 0;
theta1_dot_initial = 0;
theta2_dot_initial = 0;
theta1_dot_final = 0;
theta2_dot_final = 0;
time_initial = 0;
time_final = 10;

THETA1 = [theta1_initial theta1_dot_initial theta1_final theta1_dot_final]';
THETA2 = [theta2_initial theta2_dot_initial theta2_final theta2_dot_final]';

TIME_Coeff = [1 time_initial time_initial^2 time_initial^3; 
                     0 1 2*time_initial 3*(time_initial^2);
              1 time_final time_final^2 time_final^3; 
                     0 1 2*time_final 3*(time_final^2)];

Traj_Coeff1 = TIME_Coeff\THETA1;
Traj_Coeff2 = TIME_Coeff\THETA2;

Traj_Theta_1 = Traj_Coeff1(1) + Traj_Coeff1(2)*t + Traj_Coeff1(3)*(t^2) + Traj_Coeff1(4)*(t^3);
Traj_Theta_2 = Traj_Coeff2(1) + Traj_Coeff2(2)*t + Traj_Coeff2(3)*(t^2) + Traj_Coeff2(4)*(t^3);

fprintf("********** Desired Trajectory of Theta1 **********\n")
disp(Traj_Theta_1)

fprintf("********** Desired Trajectory of Theta2 **********\n")
disp(Traj_Theta_2)

fprintf("********** Desired Trajectory of Theta1_dot **********\n")
Traj_Theta_1_dot = diff(Traj_Theta_1);
disp(Traj_Theta_1_dot)


fprintf("********** Desired Trajectory of Theta2_dot **********\n")
Traj_Theta_2_dot = diff(Traj_Theta_2);
disp(Traj_Theta_2_dot)

fprintf("********** Desired Trajectory of Theta1_ddot **********\n")
Traj_Theta_1_ddot = diff(Traj_Theta_1,2);
disp(Traj_Theta_1_ddot)

fprintf("********** Desired Trajectory of Theta2_ddot **********\n")
Traj_Theta_2_ddot = diff(Traj_Theta_2,2);
disp(Traj_Theta_2_ddot)

