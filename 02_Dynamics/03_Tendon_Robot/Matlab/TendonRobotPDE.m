function [tip_pose_traj, tip_rot_traj] = TendonRobotPDE(motor_traj_1, motor_traj_2, T, dt)

    % Independent Parameters
    L = 0.08;          %Length (before strain)
    N = 200;          %Spatial resolution
    E = 207e9;        %Young's modulus
    shear_modulus = E/( 2*(1+0.3) );       %Shear modulus
    radiu = 0.001;   %Cross-section radius
    total_mass = 0.034;
    
    g = [-9.81;0;0];  %Gravity
    Bse = zeros(3);   %Material damping coefficients - shear and extension
    Bbt =diag([5e-4, 5e-4, 5e-4]);  %Material damping coefficients - bending and torsion
    C = diag([1e-4, 1e-4, 1e-4]);    %Square-law-drag damping coefficients
    
    %base_to_motor = 0.0518;  %length of tendon from the baseplate to the motor
    motor_shaft_radius = 0.04; % in meters
    num_tendons = 4;
    tau = [0; 0; 0; 0];
    compliance=[1e-3; 1e-3; 1e-3; 1e-3];
    num_disks=9;
    tendon_offset = 0.006;
    rot = @(theta) tendon_offset*[cos(theta); sin(theta); 0];
    r = {rot(0), rot(pi/2), rot(pi), rot(pi*3/2)}; %Tendons 90 deg apart
    
    %% Dependent parameter calculations
    area = pi*radiu^2;                             %Cross-sectional area
    J = diag([pi*radiu^4/4  pi*radiu^4/4  pi*radiu^4/2]);   %Inertia
    
    Kse = diag([shear_modulus*area, shear_modulus*area, E*area]);                %Stiffness matrix - shear and extension
    Kbt = diag([E*J(1,1), E*J(2,2), shear_modulus*J(3,3)]); %Stiffness matrix - bending and torsion
    rho = total_mass/(L*area);       %Density
    rhoA=rho*area;
    rhoAg=rho*area*g;
    %% Boundary Conditions for every time step
    p0 = [0;0;0];
    R0 = eye(3);
    q0 = [0;0;0];
    w0 = [0;0;0];
        
    %% simulation parameters
    alpha = 0;  %use BDF2
    t=0;%global time
    STEPS=T/dt; %Number of timesteps to completion
    tSteps = 0:dt:T-dt;
    
    ds = L/(N-1);                               %Grid distance (before strain)
    c0 = (1.5 + alpha) / ( dt*(1+alpha) );      %BDF-alpha coefficients
    c1 = -2/dt;
    c2 = (0.5 + alpha) / ( dt*(1+alpha) );
    d1 = alpha / (1+alpha);
    
    %% initialization
    % y:=[p;R;v;u;q;w,pib_norm] and z:=[v;u;q;w;vs;us]
    Y = [zeros(2,N); linspace(0,L,N); zeros(25,N)];
    Z = [zeros(2,N); ones(1,N); zeros(15,N)];
    for i=1:N
        Y(4:12,i)=reshape(R0,9,1);
    end
    Y_prev=Y;
    Z_prev=Z;
    % initilize the tendon actuation and the difference in motor position
    pi_b_norm = zeros(4,1);
    % initlize the tendon lengths
    l_star = [L; L; L; L];

    % init the pose and rot trajectory of the tip during the simulation
    tip_pose_traj = zeros(3,STEPS);
    tip_rot_traj  = zeros(3,3,STEPS);
    
    % fsolve parameters
    options = optimoptions('fsolve');
    options = optimoptions(options,'MaxFunctionEvaluations',1e6); 
    options = optimoptions(options,'MaxIterations',1e3); 
    
    %% Main Simulation
    % Tracking the length, stretch, slack and tension of the tendons
    tendon_stretch_t = zeros(4,STEPS);
    stretch = [];
    tendon_length_t  = zeros(4,STEPS);
    integrated_lengths = [];
    tendon_slack_t   = zeros(4,STEPS);
    slack = [];
    tendon_tension_t = zeros(4,STEPS);
   
    guessPre=fsolve(@staticBVP, zeros(6+num_tendons,1),options); %Solve static BVP w/ shooting method
    robotShow(motor_traj_1, motor_traj_2, t, tSteps);
    for tStep = 1 : STEPS
        %Set history terms
        Z_h = c1*Z+c2*Z_prev;
        Y_prev = Y;
        Z_prev = Z;
        %Midpoints are linearly interpolated for RK4
        Zh_int = 0.5*(Z_h(:,1:end-1) + Z_h(:,2:end));

        guessPre=fsolve(@staticBVP, guessPre,options);%Solve semi-discretized PDE w/ shooting 
        robotShow(motor_traj_1, motor_traj_2, t, tSteps); %3dplot
        disp(t);     %print current time
        t=t+dt;  %update time 

        % save the pose and the rotation of the tip of the TACR
        tip_pose_traj(:, tStep) = Y(1:3,end);
        tip_rot_traj(:,:,tStep) = reshape(Y(4:12,end),3,3);

        % tracking slack, stretch, length and tension of tendons
        tendon_stretch_t(:,tStep) = stretch;
        tendon_length_t(:,tStep)  = integrated_lengths;
        tendon_slack_t(:,tStep)   = slack;
        tendon_tension_t(:,tStep) = tau;
        tau

    
    end

    % plot the slack, stretch, length of tendons
    for i = 1:num_tendons
        figure()
        plot(tendon_length_t(i,:))
        hold on
        plot(tendon_stretch_t(i,:))
        hold on
        plot(tendon_slack_t(i,:))
        title(['tendon: ', num2str(i)]);
        legend('integrated lengths','stretch', 'slack');
    end
    % plot the tension of tendons
    figure()
    for i = 1:num_tendons
        plot(tendon_tension_t(i,:))
        hold on
    end
    title('tendon tension')
    legend('1','2','3','4');
    
    
    
    %ps, Rs, vs, us, pib_s_norm, nb, mb p,R,v,u
        function [y_s] = staticODE(s,y)
            %unpack
            R=reshape(y(4:12),3,3);
            v=y(13:15);
            u=y(16:18);
    
            %Setup tendon linear system
            a = zeros(3,1);
            b = zeros(3,1);
            A = zeros(3,3);
            G = zeros(3,3);
            H = zeros(3,3);
            pib_s_norm=y(24:24+num_tendons-1);
    
            for i = 1 : num_tendons
                pb_si = cross(u,r{i}) + v;
                pib_s_norm(i) = norm(pb_si);
                A_i = -hat(pb_si)^2 * (tau(i)/pib_s_norm(i)^3);
                G_i = -A_i * hat(r{i});
                a_i = A_i * cross(u,pb_si);
    
                a = a + a_i;
                b = b + cross(r{i}, a_i);
                A = A + A_i;
                G = G + G_i;
                H = H + hat(r{i})*G_i;
            end
    
            K = [A+ Kse,     G   ;
                G.' ,  H + Kbt];
    
            nb = Kse*(v - [0;0;1]);
            mb = Kbt*u;
    
            rhs = [-cross(u,nb) - R'*rho*area*g - a;
                -cross(u,mb) - cross(v,nb) - b];
    
            %Calculate ODE terms
            p_s = R*v;
            R_s = R*hat(u);
            vs_and_us = K\rhs;
            q_s=zeros(3,1);
            w_s=zeros(3,1);
    
            % pack & out y_s
            y_s=[p_s;reshape(R_s,9,1);vs_and_us;q_s;w_s;pib_s_norm];
        end
    
    %% dynamic ODE function
        function [y_s,z] = dynamicODE(y,z_h)
            %Unpack state vector
            R = reshape(y(4:12),3,3);
            v = y(13:15);
            u = y(16:18);
            q = y(19:21);
            w = y(22:24);
    
            v_h=z_h(1:3);
            u_h=z_h(4:6);
            q_h=z_h(7:9);
            w_h=z_h(10:12);
            v_sh=z_h(13:15);
            u_sh=z_h(16:18);
            %Setup tendon linear system
            a = zeros(3,1);
            b = zeros(3,1);
            A = zeros(3,3);
            G = zeros(3,3);
            H = zeros(3,3);
            pib_s_norm=zeros(num_tendons,1);
    
            for num = 1 : num_tendons
                pb_si = cross(u,r{num}) + v;
                pib_s_norm(num) = norm(pb_si);
                A_i = -hat(pb_si)^2 * (tau(num)/pib_s_norm(num)^3);
                G_i = -A_i * hat(r{num});
                a_i = A_i * cross(u,pb_si);
    
                a = a + a_i;
                b = b + cross(r{num}, a_i);
                A = A + A_i;
                G = G + G_i;
                H = H + hat(r{num})*G_i;
            end
    
            K = [A+ Kse+c0*Bse,     G   ;
                G.' ,  H + Kbt+c0*Bbt];
    
            v_t = c0*v + v_h;
            u_t = c0*u + u_h;
            q_t = c0*q + q_h;
            w_t = c0*w + w_h;
    
            nb = Kse*(v - [0;0;1])+Bse*v_t;
            mb = Kbt*u+Bbt*u_t;
    
            rhs = [-a + rhoA*(cross(w,q)+q_t) + C*q.*norm(q) - R'*rhoAg - cross(u,nb)-Bse*v_sh;
                -b + rho*cross(w,J*w)+rho*J*w_t - cross(v,nb)-cross(u,mb)-Bbt*u_sh];
    
            %ODEs
            p_s = R*v;
            R_s = R*hat(u);
            vs_and_us=K\rhs;
            q_s = v_t - hat(u)*q + hat(w)*v;
            w_s = u_t - hat(u)*w;
    
            %Pack state vector derivative
            z=[v;u;q;w;vs_and_us];
            y_s = [p_s; reshape(R_s,9,1); vs_and_us;q_s;w_s;pib_s_norm];
    
        end
    
    %% static BVP shoot
        function distal_error = staticBVP(guess)
            pi_b_norm = tendon_q(motor_traj_1, motor_traj_2, motor_shaft_radius, t, dt);
            n0 = guess(1:3);
            v0 = Kse\n0 + [0;0;1];
            u0 = guess(4:6);
            tau=max(guess(7:7+num_tendons-1),0); %Pull force, which is the largest when compared with 0
            slack=-min(guess(7:7+num_tendons-1),0);
            
            s=linspace(0,L,N);
            %ODE45
            %  y0=[p0; reshape(R0,9,1); v0; u0; q0; w0; pi_b_norm];
            %  Sol=ode45(@staticODE, s, y0);
            %  [Y,y_s] =deval(Sol,s); %obtain Y and derivatives y_s
            %  for j = 1 : N-1
            %    Z(1:6,j)=Y(13:18,j);
            %    Z(13:18,j)=y_s(13:18,j);
            %  end
    
            %Euler's method
            Y(:,1) = [Y(1:12,1); v0; u0; q0; w0; pi_b_norm]; %initial y0
            for j = 1 : N-1
                [y_s] = staticODE(s(j),Y(:,j));
                Y(:,j+1)=Y(:,j)+ds*y_s;
                Z(1:6,j)=Y(13:18,j);
                Z(13:18,j)=y_s(13:18);
            end
    
            %Find the internal forces in the backbone prior to the final plate
            vL = Y(13:15,end);
            uL = Y(16:18,end);
            nb_L = Kse*(vL - [0;0;1]);
            mb_L = Kbt*uL;
    
            %Find the equilibrium error at the tip, considering tendon forces
            force_error = -nb_L;
            moment_error = -mb_L;
    
            for index = 1 : num_tendons
                pb_si = cross(uL,r{index}) + vL;
                Fb_i = -tau(index)*pb_si/norm(pb_si);
                force_error = force_error + Fb_i;
                moment_error = moment_error + cross(r{index}, Fb_i);
            end
    
            %Find the length violation error
            integrated_lengths=Y(25:25+num_tendons-1,end);
            %l_star=L+base_to_motor+Z(t);
            
            stretch = l_star.*(compliance.*tau);
            length_error = integrated_lengths + slack - (l_star + stretch);
            distal_error = [force_error; moment_error;length_error];
        end
    %% dynamic BVP shoot
        function distal_error = dynamicBVP(guess)
            pi_b_norm = tendon_q(motor_traj_1, motor_traj_2, motor_shaft_radius, t, dt);
    
            n0 = guess(1:3);
            v0 = Kse\n0 + [0;0;1];
            u0 = guess(4:6);
            tau=max(guess(7:7+num_tendons-1),0); %Pull force, which is the largest when compared with 0
            slack=-min(guess(7:7+num_tendons-1),0);
    
            %Y(:,1) = [Y(1:24,1); pi_b_norm]; %initial y0
            Y(:,1) = [Y(1:12,1); v0; u0; q0; w0; pi_b_norm]; %initial y0
            
            for j = 1 : N-1
                % Euler's method
                [y_s,Z(:,j)] = dynamicODE(Y(:,j),Z_h(:,j));
                Y(:,j+1)=Y(:,j)+ds*y_s;
                %OED45
                % Yj = Y(:,j);
                % [k1,Z(:,j)]=dynamicODE(Yj,Z_h(:,j));
                % [k2,~]=dynamicODE(Yj+k1*ds/2,Zh_int(:,j));
                % [k3,~]=dynamicODE(Yj+k2*ds/2,Zh_int(:,j));
                % [k4,~]=dynamicODE(Yj+k3*ds,Z_h(:,j+1));
                % Y(:,j+1) = Yj + ds*(k1 + 2*(k2+k3) + k4)/6;
            end
            
            %Find the internal forces in the backbone prior to the final plate
            vL = Y(13:15,end);
            uL = Y(16:18,end);
            vL_t=c0*vL + Z_h(1:3,end);
            uL_t=c0*vL + Z_h(4:6,end);
    
            nb_L = Kse*(vL - [0;0;1])+Bse*vL_t;
            mb_L = Kbt*uL + Bbt*uL_t;
    
            %Find the equilibrium error at the tip, considering tendon forces
            force_error = -nb_L;
            moment_error = -mb_L;
    
            for index = 1 : num_tendons
                pb_si = cross(uL,r{index}) + vL;
                Fb_i = -tau(index)*pb_si/norm(pb_si);
                force_error = force_error + Fb_i;
                moment_error = moment_error + cross(r{index}, Fb_i);
            end
            %Find the length violation error
            integrated_lengths=Y(25:25+num_tendons-1,end);
            %l_star=L+base_to_motor+Z(t);
            stretch = l_star.*(compliance.*tau);
            length_error = integrated_lengths + slack - (l_star + stretch);
            distal_error = [force_error; moment_error;length_error];
    
        end
    
    
    %% other func
        function skew_symmetric_matrix = hat(y)
            skew_symmetric_matrix = [  0   -y(3)  y(2) ;
                                    y(3)   0   -y(1) ;
                                    -y(2)  y(1)   0  ];
        end
    
        function robotShow(motor_traj_1, motor_traj_2, t, tSteps)
            centerline=Y(1:3,:);
            tendonlines = zeros(3*num_tendons, size(centerline,2));
            for i = 1 : size(centerline,2)
                p_show = Y(1:3,i);
                R_show = reshape(Y(4:12,i),3,3);
                for j = 1 : num_tendons
                    tendonlines(3*j-2 : 3*j, i) = p_show + R_show*r{j};
                end
            end
            disks = zeros(3,4*num_disks);
            for i = 1 : num_disks
                j = round(size(centerline,2) * i / num_disks);
                p_show = Y(1:3,j);
                R_show = reshape(Y(4:12,j),3,3);
                disks(1:3, 4*i-3:4*i-1) = R_show;
                disks(1:3, 4*i) = p_show;
            end
            Visualize(centerline,tendonlines,num_tendons, disks,num_disks, motor_traj_1, motor_traj_2, t, tSteps);
        end
end