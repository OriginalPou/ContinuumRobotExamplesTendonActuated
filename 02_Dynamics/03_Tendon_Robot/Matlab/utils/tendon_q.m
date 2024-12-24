%       Function tendon_q:
%           Computes the actuation of each tendon, where tendons are only pulled on released
%           Like it is the case for Tendon Actuated Continuum Robots         
%
%           Parameters
%           ----------
%                           motor_traj -> the desired trajectory of the two 
%                                 controlling each of the two pairs
%                                 of antagonistic cables (in radians)
%                           motor_shaft_radius -> the radius of the motor shaft
%                                  tendon_q = motor_traj(rad) * motor_shaft_radius
%                           t -> the current time of the simulation
%                           dt -> time step
%           
%           Returns
%           -------
%                           tendon_q -> a (4,1) array of the actuation of each tendon
%
%           Note
%           ----
%                           
%           - This function allows the tendon actuated continuum robot to always
%           be actuated by pulling on one of the antagonistic cables
%
%           - The tendons are pulled when the motor turns in one direction,
%           When the motor changes direction, the tendon is released until
%           goes back to its initial configuration and then its antagonistic
%           tendon is pulled
%
%           - When the motor position is positive, the tendons 1 and 2 are pulled
%           otherwise, the tendons 3 and 4 are pulled

function tendon_q = tendon_q(motor_1_traj, motor_2_traj, motor_shaft_radius, t, dt)
    
    % compute the index corresponding to the current time
    idx = max(1, round(t/dt));
    
    % initilize the tendon actuation
    tendon_q = zeros(4,1);

    % the first pair of antagonistic tendons 1-3
    if motor_1_traj(idx)>=0
        tendon_q(1) = motor_1_traj(idx);
    else
        tendon_q(3) = -motor_1_traj(idx);
    end

    % the second pair of antagonistic tendons 2-4
    if motor_2_traj(idx)>=0
        tendon_q(2) = motor_2_traj(idx);
    else
        tendon_q(4) = -motor_2_traj(idx);
    end 

    tendon_q = tendon_q * motor_shaft_radius; % in meters
end