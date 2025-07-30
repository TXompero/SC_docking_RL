% RL for approach maneuver in SC docking, simplified Environment
clear; clc; close all

%% Set parameters

% Orbital and chief parameters
GM_e = 398600;     % [km^3/s^2]

a_chief = 6780;     % [km]
e_chief = 0;        % [-]
i_chief = 0;        % [deg]
OM_chief = 0;       % [deg]
om_chief = 0;       % [deg]
th_chief = 0;       % [deg]
kep_chief = [a_chief, e_chief, i_chief, OM_chief, om_chief, th_chief];

n_chief = sqrt(GM_e / a_chief^3);

% Control parameters
Ts = 1;
Tf = 3600;
max_acc = 0.05;

% Define step function weights and starting parameters for the reset function
y_start = 100;              % [m]
y_target = 10;              % [m]
r_start = 50;               % [m]
r_target = 5;               % [m]
cone_angle = 30*pi/180;     % [rad]
min_velocity_target = 0.01; % [m/s]

init_rand_vel = 0.1;        % [m/s]
init_vel_range = 0.01;      % [m/s]
vel_angle = 45*pi/180;      % [rad]

%% Create the environment

% Define observation and actuator informations
obsInfo = rlNumericSpec([6 1]);
obsInfo.Name = "Cartesian States";
obsInfo.Description = 'x, y, z, vx, vy, vz';

actInfo = rlNumericSpec([3 1]);
actInfo.Name = "Impulsive Delta v Action";
actInfo.LowerLimit = -max_acc;
actInfo.UpperLimit = max_acc;

% Set integration options for the step function
int_opts = odeset('AbsTol',1e-12,'RelTol',1e-12);

% Define environment for the training 
env = rlFunctionEnv(obsInfo,actInfo,...
    @(Action,State) StepFunction(Action, State, Ts, n_chief, max_acc, cone_angle, y_start,...
    r_target, y_target, min_velocity_target, int_opts),...
    @() ResetFunction(y_start,r_start,vel_angle,init_rand_vel,init_vel_range));

numObs = obsInfo.Dimension(1);
numAct = actInfo.Dimension(1);

%% Create neural networks for Actor, Critic

% Initialize NN
initOpts = rlAgentInitializationOptions(NumHiddenUnit=128);

% Define settings
criticOptions = rlOptimizerOptions( ...
    Optimizer="adam", ...
    LearnRate=0.001, ...
    GradientThreshold=1);

actorOptions = rlOptimizerOptions( ...
    Optimizer="adam", ...
    LearnRate=0.0005, ...
    GradientThreshold=1);

% Additional optinal settings
agentOptions = rlTD3AgentOptions;
agentOptions.DiscountFactor = 0.999;
agentOptions.TargetSmoothFactor = 1e-3;
agentOptions.MiniBatchSize = 128;
agentOptions.CriticOptimizerOptions = criticOptions;
agentOptions.ActorOptimizerOptions = actorOptions;

% Exploration noise settings
agentOptions.ExplorationModel.StandardDeviationMin = 0;
agentOptions.ExplorationModel.StandardDeviation = 0.4;
agentOptions.ExplorationModel.StandardDeviationDecayRate = 0.0005;

% Fix randomness (optional), uncomment line below if wanted
rng(0,"twister");

% Create agent
agent = rlTD3Agent(obsInfo,actInfo,initOpts,agentOptions);

%% Train agent

% Training options
% StopTrainingCriteria: "EpisodeCount"       : Stop at the n episode evaluated
%                       "AverageReward"      : Stop at the specified average reward achieved
%                       "EvaluationStatistic": Stop when the evaluator reward value is the specified one
maxepisodes = 5000;
maxsteps = ceil(Tf/Ts);
trainingOptions = rlTrainingOptions(...
    MaxEpisodes=maxepisodes,...
    MaxStepsPerEpisode=maxsteps,...
    ScoreAveragingWindowLength=20,...
    Verbose=false,...
    Plots="training-progress",...
    StopTrainingCriteria="AverageReward",...
    StopTrainingValue=150, ...
    UseParallel=true);

% Agent evaluator (Optional), used to obtain episodes with no exploration noise
% EvaluationStatisticType: "MinEpisodeReward" : Of all evaluations, consider the minimum one
%                          "MeanEpisodeReward": Consider the mean value of all evaluator episode reward
evl = rlEvaluator(NumEpisodes=8,EvaluationFrequency=250,EvaluationStatisticType="MinEpisodeReward");

% Perform train
trainingStats = train(agent,env,trainingOptions,Evaluator=evl);

%% Evaluate agent solution

verify_agent;
