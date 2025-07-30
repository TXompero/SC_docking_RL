%% Build enviroment plot

init_r = 50;
init_y = 100;

clear length

th = linspace(0, 2*pi, 1e3);

% Plot the sphere
figure(1)
plot3(0,0,0,'ko','LineWidth',2,'DisplayName','Target spacecraft')
hold on
axis equal
xlabel('LVLH x axis [m]'); ylabel('LVLH y axis [m]'); zlabel('LVLH z axis [m]');
fill3(r_start*cos(th),-y_start*ones(length(th)),r_start*sin(th),'r','FaceAlpha',0.2,'DisplayName','Initial position region')
plot3(r_target*cos(th),-y_target*ones(length(th)),r_target*sin(th),'k','LineWidth',1.5,'DisplayName','Target circle position')
title('Trajectory visualization')
grid on

% Plot cone
height = init_y+5;        
radius = height*cosd(30);         
n = 100;                      % Resolution

% Generate angle and height arrays
h = linspace(0, height, n);
[Theta, H] = meshgrid(th, h);

% Radius tapers linearly from tip to base
R = (H / height) * radius;

% Cone in local coordinates (Y-axis as height direction)
X = R .* cos(Theta);
Z = R .* sin(Theta);
Y = -H;    % Negative to point along -Y direction

% Plot
surf(X, Y, Z, 'FaceColor', 'k', 'EdgeColor', 'none','FaceAlpha',0.1)

%% Evaluate episode trajectories

% Set number of episodes
it = 100;

% Initialize parameters
length = 100;
xx0_matr = zeros(6,length);
flg_matr = zeros(2,length);
fail = zeros(it,1);
reward = zeros(it,1);

vel_angle = zeros(it,1);
vel_treshold = zeros(it,1);
vel_init = zeros(it,1);
final_pos = zeros(it,2);

clear length

% Set propagation parameters
% GM_e = 398600;     % [km^3/s^2]
% a_chief = 6780;    % [km]
% n_chief = sqrt(GM_e / a_chief^3);
t_end = 0;

figure(2)
subplot(3,1,1); grid on; hold on
title('Spacecraft velocity x axis')
xlabel('Time [s]')
ylabel('Velocity [m/s]')
subplot(3,1,2); grid on; hold on
title('Spacecraft velocity y axis')
xlabel('Time [s]')
ylabel('Velocity [m/s]')
subplot(3,1,3); grid on; hold on
title('Spacecraft velocity z axis')
xlabel('Time [s]')
ylabel('Velocity [m/s]')

figure(3)
title('Accelerations')
subplot(3,1,1); grid on; hold on
title('Spacecraft acceleration x axis')
xlabel('Time [s]')
ylabel('Accel. [m/s^2]')
subplot(3,1,2); grid on; hold on
title('Spacecraft acceleration y axis')
xlabel('Time [s]')
ylabel('Accel. [m/s^2]')
subplot(3,1,3); grid on; hold on
title('Spacecraft acceleration z axis')
xlabel('Time [s]')
ylabel('Accel. [m/s^2]')

for i = 1:it

    % Evaluate start condition
    xx0 = reset(env);
    
    % Initialize saving vector
    xx_all = [];

    if it == 1
        xx0_rand = xx0;
    end

    % Initialize additional parameters
    t_end = 0;
    isDone = false;
    
    while ~isDone && t_end < Tf
        
        [current_action, current_info] = getAction(agent,xx0);
        [tt_prop,xx_prop] = ode113(@(t,x) CW_eom(x,t,current_action{1},n_chief),[0 Ts],xx0,int_opts);
    
        xx_all = [xx_all; xx_prop];

        % Step in environment
        [nextObs, r, isDone, nextState] = step(env, current_action{1});

        % Cumulative reward
        reward(i) = reward(i) + r;

        % Save initial velocity
        if t_end == 0
            vel_init(i) = norm(xx0(4:6));
        end

        % For first iteration plot velocity and acceleration evolution over time
        if i == 1
            
            figure(2)
            subplot(3,1,1); plot(tt_prop+t_end,xx_prop(:,4),'b','LineWidth',1.5)
            subplot(3,1,2); plot(tt_prop+t_end,xx_prop(:,5),'b','LineWidth',1.5)
            subplot(3,1,3); plot(tt_prop+t_end,xx_prop(:,6),'b','LineWidth',1.5)

            figure(3)
            subplot(3,1,1); plot(t_end+[0, Ts],current_action{1}(1)*[1, 1],'b','LineWidth',1.5)
            subplot(3,1,2); plot(t_end+[0, Ts],current_action{1}(2)*[1, 1],'b','LineWidth',1.5)
            subplot(3,1,3); plot(t_end+[0, Ts],current_action{1}(3)*[1, 1],'b','LineWidth',1.5)

        end
    
        xx0 = (xx_prop(end,:))';
        t_end = t_end + tt_prop(end);

        % Verify consistency between enviroment and propagation
        if max(abs(xx0 - nextState)) > 1e-3
            error('Inconsistency')
        end
    
    end
    
    % Save fail conditions
    if sqrt(xx0(1)^2 + xx0(3)^2) > r_target
        fail(i) = 1;
    end

    if abs(xx0(2)) > y_target
        fail(i) = 3;
    end

    % Plot trajectory
    figure(1)
    plot3(xx_all(:,1),xx_all(:,2),xx_all(:,3),'b','LineWidth',1.5)

    % Save parameters for additional plots
    vel_fin = xx0(4:6);
    vel_angle(i) = acos(dot(vel_fin/norm(vel_fin),[0; 1; 0]));
    vel_treshold(i) = norm(vel_fin);
    final_pos(i,:) = [xx0(1), xx0(3)];

end

% Fill additional plot
figure(4)
grid on; hold on
plot(1:it,vel_treshold,'LineWidth',1.5)
xlabel('Iteration')
ylabel('Final velocity [m/s]')
title('Final velocity')

figure(5)
grid on; hold on
plot(1:it,vel_init,'LineWidth',1.5)
xlabel('Iteration')
ylabel('Initial velocity [m/s]')

figure(6)
grid on; hold on
plot(1:it,reward,'LineWidth',1.5)
xlabel('Iteration')
ylabel('Episode reward')

figure(7)
grid on; hold on
plot(r_target*cos(th),r_target*sin(th),'k--','LineWidth',1.5)
legend('r = 5 m','Autoupdate','off')
plot(final_pos(:,1),final_pos(:,2),'*','LineWidth',1.5)
xlabel('LVLH x axis [m]')
ylabel('LVLH z axis [m]')
title('Final position, y = -10 m')
axis equal