function [InitialObservation, InitialState] = ResetFunction(y_init, rad_init, v_angle, vel_midpoint, vel_range)

rand_numbers = rand(5,1);

% Generate random position
rad = rad_init*rand_numbers(1);   
angle1 = 2*pi*rand_numbers(2);

rr = [rad*cos(angle1); -y_init; rad*sin(angle1)];

% Generate random velocity and rotate
vel = 2*vel_range*rand_numbers(3) + vel_midpoint - vel_range;   % m/s
angle1 = (v_angle*2)*rand_numbers(4)-v_angle; angle2 = (v_angle*2)*rand_numbers(5) - v_angle;
angle1 = angle1*pi/180; angle2 = angle2*pi/180;

vv = [vel*sin(angle1)*cos(angle2); vel*cos(angle1)*cos(angle2); vel*sin(angle2)];

% Generate output
InitialState = [rr; vv];
InitialObservation = InitialState;