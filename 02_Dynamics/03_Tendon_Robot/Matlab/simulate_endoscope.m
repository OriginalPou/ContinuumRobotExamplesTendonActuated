close all
clear all
clc
%%
% the trajectory of the motor controlling the tendon 1-3

function motor_traj_1 = motor_traj_1(t)
    phase_shift = 0;
    speed = [0.8, 1.6];
    if t < phase_shift
        motor_traj_1 = 0;
    elseif t >= phase_shift & t < phase_shift + (2 * pi)/ speed(1)
        motor_traj_1 = sin(speed(1) * (t-phase_shift)); % in rad
    else
        phase_shift_ = phase_shift + (2 * pi)/ speed(1);
        motor_traj_1 = sin(speed(2) * (t-phase_shift_)); % in rad
    end    
end
phase_shift_act_1 = 0;
speed_act_1 = [0.8, 1.6];
%%
% the trajectory of the motor controlling the tendon 2-4
function motor_traj_2 = motor_traj_2(t)
    phase_shift = pi/5;
    speed = [1.5, 0.7];
    if t < phase_shift
        motor_traj_2 = 0;
    elseif t >= phase_shift & t < phase_shift + (2 * pi)/ speed(1)
        motor_traj_2 = sin(speed(1) * (t-phase_shift)); % in rad
    else
        phase_shift_ = phase_shift + (2 * pi)/ speed(1);
        motor_traj_2 = sin(speed(2) * (t-phase_shift_)); % in rad
    end 
end
phase_shift_act_2 = pi/5;
speed_act_2 = [1.5, 0.7];
%%
T = 20; % simulation duration in seconds
dt = 0.05; % timestep in seconds
t_span  = 0:dt:T-dt;
STEPS=T/dt; %Number of timesteps to completion
%%
% Define the trajectory at the proximal side 
prox_motor_traj_1 = @motor_traj_1; 
prox_motor_traj_2 = @motor_traj_2;
%% Simulation without bending shaft
% Ignore the flexible shaft
dist_motor_traj_1 = zeros(STEPS,1);
dist_motor_traj_2 = zeros(STEPS,1);
i = 1;
for time = 0:dt:T-dt
    dist_motor_traj_1(i) = prox_motor_traj_1(time);
    dist_motor_traj_2(i) = prox_motor_traj_2(time);
    i = i+1;
end
%%
[tip_pos, tip_rot] = TendonRobotPDE(dist_motor_traj_1, dist_motor_traj_2, T, dt);

%% Simulation with hysteresis from the bending shaft of the endoscope
% Coleman Hodgon model params
c = 1;
beta = -0.91;
alpha = 10.1;
% compute the motor trajectory at the distal end
dist_motor_traj_1_hyst = coleman_hodgon_model(dist_motor_traj_1, t_span, alpha, beta, c); % in rad
dist_motor_traj_2_hyst = coleman_hodgon_model(dist_motor_traj_2, t_span, alpha, beta, c); % in rad
%% plot the motor trajectory at the distal end
figure()
plot(t_span, dist_motor_traj_1_hyst);
hold on
plot(t_span, dist_motor_traj_1);
hold on
plot(t_span, dist_motor_traj_2_hyst);
hold on
plot(t_span, dist_motor_traj_2);
legend(["act 1 hyst","act 1", "act 2 hyst", "act 2"])
%%
[tip_pos_hyst, tip_rot_hyst] = TendonRobotPDE(dist_motor_traj_1_hyst, dist_motor_traj_2_hyst, T, dt);
%% Plotting results
figure()
%plot(tip_pos(1,:))
hold on
%plot(tip_pos(2,:))
hold on
%plot(tip_pos(3,:))
hold on
plot(tip_pos_hyst(1,:))
hold on
plot(tip_pos_hyst(2,:))
hold on
plot(tip_pos_hyst(3,:))
hold on
title('tip position')
xlabel('time')
legend('x','y','z', 'x-hyst','y-hyst','z-hyst');
%%
figure()
plot(t_span, 0.04*dist_motor_traj_1_hyst);
hold on
plot(t_span, 0.04*dist_motor_traj_1);
hold on
plot(t_span, 0.04*dist_motor_traj_2_hyst);
hold on
plot(t_span, 0.04*dist_motor_traj_2);
legend(["act 1 hyst","act 1", "act 2 hyst", "act 2"])
%% export data
% save('data/pf_tests/tip_pos_hyst.mat','tip_pos_hyst');
% save('data/pf_tests/tip_rot_hyst.mat','tip_rot_hyst');
% %% export data without bending part
% save('data/frames/tip_pos.mat', 'tip_pos')
% save('data/frames/tip_rot.mat', 'tip_rot')
% %% export motor traj
% q1 = dist_motor_traj_1_hyst;
% q2 = dist_motor_traj_2_hyst;
% save('data/pf_tests/q1_hyst.mat', 'q1')
% save('data/pf_tests/q2_hyst.mat', 'q2')
% %%
% q1 = dist_motor_traj_1;
% q2 = dist_motor_traj_2;
% save('data/pf_tests/q1.mat', 'q1')
% save('data/pf_tests/q2.mat', 'q2')

