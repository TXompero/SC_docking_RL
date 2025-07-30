close all; clearvars; clc

addpath('Auxiliary\')

%% Define enviroment and parameters

% Parameter of planet, chief and chaser

GM_e = 398600;     % [km^3/s^2]

a_chief = 6780;     % [km]
e_chief = 0;        % [-]
i_chief = 30;        % [deg]
OM_chief = 119.3;       % [deg]
om_chief = 110.4;       % [deg]
th_chief = 160;       % [deg]
kep_chief = [a_chief, e_chief, i_chief, OM_chief, om_chief, th_chief];

rand_number = rand(5,1);
a_chaser = a_chief + 0.4*rand_number(1) - 0.2;
e_chaser = 0;
i_chaser = i_chief + 0.04*rand_number(2) - 0.02;
OM_chaser = OM_chief + 0.04*rand_number(3) - 0.02;
om_chaser = om_chief + 0.04*rand_number(4) - 0.02;
th_chaser = th_chief + 0.04*rand_number(5) - 0.02;
kep_chaser = [a_chaser, 0, i_chaser, OM_chaser, om_chaser, th_chaser];

% Compute additional parameters
n_chief = sqrt(GM_e / a_chief^3);

R_first = [cosd(OM_chief), sind(OM_chief), 0; -sind(OM_chief), cosd(OM_chief), 0; 0, 0, 1];
R_second = [1, 0, 0; 0, cosd(i_chief), sind(i_chief); 0, -sind(i_chief), cosd(i_chief)];
R_third = [cosd(om_chief+th_chief), sind(om_chief+th_chief), 0; -sind(om_chief+th_chief), cosd(om_chief+th_chief), 0; 0, 0, 1];
R_eci2lvlh = (R_first*R_second*R_third)';


%% Verify parameters, obtain positions

xx_0_chief = kep2car(kep_chief, GM_e);
xx_0_chaser = kep2car(kep_chaser, GM_e);

ode_opts = odeset('AbsTol',1e-10,'RelTol',1e-10);
[tt_chief,xx_chief] = ode45(@(t,x) tbp_eom(x,t,GM_e),[0 2*pi/n_chief],xx_0_chief, ode_opts);

LVLH_x = a_chief/4*xx_0_chief(1:3) / norm(xx_0_chief(1:3));
LVLH_z = a_chief/4*cross(xx_0_chief(1:3),xx_0_chief(4:6))/norm(cross(xx_0_chief(1:3),xx_0_chief(4:6)));
LVLH_y = a_chief/4*cross(LVLH_z,LVLH_x) / norm(cross(LVLH_z,LVLH_x));

% Build figure
Terra_3D
plot3(xx_chief(:,1),xx_chief(:,2),xx_chief(:,3),'LineWidth',1.5)
hold on; grid on; 
% Plot LVLH ref. frame
plot3(xx_chief(1,1),xx_chief(1,2),xx_chief(1,3),'ro')
quiver3(xx_chief(1,1),xx_chief(1,2),xx_chief(1,3),LVLH_x(1),LVLH_x(2), LVLH_x(3),'r','LineWidth',1.5,'MaxHeadSize',2);
quiver3(xx_chief(1,1),xx_chief(1,2),xx_chief(1,3),LVLH_y(1),LVLH_y(2), LVLH_y(3),'c','LineWidth',1.5,'MaxHeadSize',2);
quiver3(xx_chief(1,1),xx_chief(1,2),xx_chief(1,3),LVLH_z(1),LVLH_z(2), LVLH_z(3),'g','LineWidth',1.5,'MaxHeadSize',2);
quiver3(0,0,0,xx_0_chief(1),xx_0_chief(2),xx_0_chief(3),1,'k--','LineWidth',1.5)
plot3(xx_0_chaser(1),xx_0_chaser(2),xx_0_chaser(3),'k*','LineWidth',1.5)
axis equal
xlabel('ECI x-axis [km]')
ylabel('ECI y-axis [km]')
zlabel('ECI z-axis [km]')
legend('','Target orbit','','LVLH x-axis','LVLH y-axis','LVLH z-axis')

% Get difference vectors
rr_chaser_lvlh = R_eci2lvlh*(xx_0_chaser(1:3) - xx_0_chief(1:3));
vv_chaser_lvlh = R_eci2lvlh*(xx_0_chaser(4:6) - xx_0_chief(4:6)) + n_chief*[rr_chaser_lvlh(2); rr_chaser_lvlh(1); 0];

xx0_chaser_lvlh = [rr_chaser_lvlh; vv_chaser_lvlh];

%% Verify CW equations

PHI_t_f = get_STD_CW(n_chief,5e2);
xx_final = PHI_t_f*xx0_chaser_lvlh;

PHI_rr = PHI_t_f(1:3,1:3);
PHI_rv = PHI_t_f(1:3,4:6);
delta_v0 = -(PHI_rr*xx0_chaser_lvlh(1:3) + PHI_rv*xx0_chaser_lvlh(4:6));
delta_v0 = PHI_rv \ delta_v0;

[tt_chaser_lvlh,xx_chaser_lvlh] = ode45(@(t,x) CW_eom(x,t,zeros(3,1),n_chief),[0 5e2],...
    xx0_chaser_lvlh+[zeros(3,1); delta_v0], ode_opts);

figure(2)
plot3(xx_chaser_lvlh(:,1),xx_chaser_lvlh(:,2),xx_chaser_lvlh(:,3),'LineWidth',1.5)
hold on; grid on
plot3(0,0,0,'ko','LineWidth',1.5)
plot3(xx_chaser_lvlh(1,1),xx_chaser_lvlh(1,2),xx_chaser_lvlh(1,3),'ro','LineWidth',1.5)
legend('LVLH chaser trajectory','Chief','Start')

%% Visualize Approach zone and starting conditions

% Cone parameters
height = 100;        % Cone height
radius = height*cosd(30);         % Cone base radius
n = 100;             % Resolution

% Generate angle and height arrays
theta = linspace(0, 2*pi, n);
h = linspace(0, height, n);
[Theta, H] = meshgrid(theta, h);

% Radius tapers linearly from 0 (tip) to base
R = (H / height) * radius;

% Cone in local coordinates (Y-axis as height direction)
X = R .* cos(Theta);
Z = R .* sin(Theta);
Y = -H;    % Negative to point along -Y direction

% Plot
figure(3)
surf(X, Y, Z, 'FaceColor', 'k', 'EdgeColor', 'none','FaceAlpha',0.1)
hold on; grid on
plot3(0,0,0,'ko','LineWidth',2)
plot3(5*cos(theta),-10*ones(length(theta)),5*sin(theta),'k','LineWidth',1.5)
legend('Approach cone','Target S/C','Target position','AutoUpdate','off')
fill3(50*cos(theta),-100*ones(length(theta)),50*sin(theta),'r','FaceAlpha',0.5)
%legend('Initial position circle','AutoUpdate','off')
axis equal
xlabel('LVLH x-axis [m]')
ylabel('LVLH y-axis [m]')
zlabel('LVLH z-axis [m]')
title('Actual enviroment')



