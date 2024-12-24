clear all
close all


function displace=Z_t(t)
    zA = 0;          %  z(m)  Ramp Motor Input:
    zB = -0.025;     % zA ^ ____           ______
    t1 = 1;          %    |     \         /
    t2 = 2;          % zB |      \_______/
    t3 = 3;         
 %    -------------------------> t(s)
    t4 = 4;          %         t1 t2    t3 t4
    if t>t1 && t<=t2
        displace=zA + (zB-zA)*(t-t1)/(t2-t1); %Ramp lower
    elseif t>t2 && t<=t3
        displace=zB;                         %Low state
    elseif t>t3 && t<=t4
        displace=zB + (zA-zB)*(t-t3)/(t4-t3); %Ramp higher
    else
        displace=zA;  %High state
    end
end

dt = 0.01;
time = -2*pi:dt:2*pi;
motor_traj = [];

n_steps = length(time);

%%
% 
% motor_q_one_direction = [0];
% for i = 2:n_steps
%     motor_traj_dt = motor_traj(i)-motor_traj(i-1);
%     if (motor_q_one_direction(end) == 0) | (sign(motor_traj_dt) == sign(motor_q_one_direction(end)))
%         motor_q_one_direction = [motor_q_one_direction, motor_q_one_direction(end)+motor_traj_dt];
%     else
%         motor_q_one_direction = [motor_q_one_direction,0];
%     end
% end
% 
% figure();
% plot(t,motor_q_one_direction);
% hold on
% plot(t, motor_traj);

%%
tendon_q = zeros(4,1);
motor_dq = zeros(2,1);

flags = [];
tendon_qs = [];
compliances = [];

for t = -2*pi:dt:2*pi
    [tendon_q, compliance, motor_dq, direction_flag] = tendon_ex_length(@Z_t, @Z_t, tendon_q, motor_dq, t, dt);
    
    flags = [flags, direction_flag];
    tendon_qs = [tendon_qs, tendon_q];
    compliances = [compliances, compliance];

    motor_traj = [motor_traj, Z_t(t)];


end
%%
figure();
for i = [2,4]
    plot(tendon_qs(i,:));
    hold on
end
%legend('1','2','3','4')
hold on
plot(flags(1,:));
hold on
plot(motor_traj);
