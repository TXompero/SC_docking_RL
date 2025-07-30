function [NextObs,Reward,IsDone,NextState] = StepFunction(Action,State,Ts,n_chief,...
    MaxImpulse,coneAngle,y_max,r_circ_target,y_target,min_velocity,odeOpts)

% Retrieve state
rr = State(1:3);
vv = State(4:6);

% Bound action, safe step
Action(Action > MaxImpulse) = MaxImpulse;
Action(Action < -MaxImpulse) = -MaxImpulse;

% Compute next step
[~,xx_prop] = ode113(@(t,x) CW_eom(x,t,Action(:),n_chief),[0 Ts],[rr; vv],odeOpts);

NextState = xx_prop(end,1:6);
NextState = NextState';
NextObs = NextState;

% Compute reward
err_pos_d = (abs(NextState(2)))/(y_max);
err_pos_circ = norm([NextState(1); 0; NextState(3)])/r_circ_target;

Reward = (1 - err_pos_d) + (1 - err_pos_circ);

% STop if SC has reached target y coordinate
IsDone = abs(NextState(2)) < y_target;

% Additional reward to incentivize end inside target circle
if IsDone && norm([NextState(1); 0; NextState(3)]) < r_circ_target
    Reward = Reward + 3; 
end

% Verify other stop conditions (approach cone, minimum velocity, SC aument distance with target, outside taget circle)
vec1 = NextState(1:3) / norm(NextState(1:3));
vec2 = [0; -1; 0];
angle = acos(dot(vec1,vec2));

if angle > coneAngle || NextState(2) < State(2) || norm(NextState(4:6)) < min_velocity ...
        || (IsDone && norm([NextState(1); 0; NextState(3)]) > r_circ_target) 
    Reward = Reward - 100;
    IsDone = true;
end
