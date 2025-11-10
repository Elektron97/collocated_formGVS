%%% Collocated Form Control %%%
clear all
close all
clc

%% Load and Startup SoRoSim
% Clean StartUp
diff_sorosim_path = fullfile("SoRoSim", "Differentiable_SoRoSim");
cd(diff_sorosim_path)
startup

% Switch again to the current directory
[current_path, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(current_path)

%% Load Data
robot_name = "rsip";
% robot_name = "conical_hsupport";

% mat ext
file_name = "robot_linkage";
mat_ext = ".mat";

% Load Robot and Data
load(fullfile("robots", robot_name, "robot_linkage" + mat_ext));

%% Update Robot
% Camera Position
T1.PlotParameters.Light = false;
T1.PlotParameters.ClosePrevious = false;

% Axes Limits
T1.PlotParameters.XLim = [-T1.VLinks.L, T1.VLinks.L];
T1.PlotParameters.YLim = [-T1.VLinks.L, T1.VLinks.L];

% Colors
blue_sofft = "#086788";
red_target = "#f06543";
grey_mid = "#858583";
T1.VLinks.color = hex2rgb(blue_sofft);

% Camera Position
T1.PlotParameters.CameraPosition = [-0.0316   -0.0001   -6.9004];

% CC Segment
T1.CVRods{1}(2).Phi_odr = zeros(6, 1);

% Update Linkage
T1 = T1.Update();

% Damping Joint
if T1.CVRods{1}(1).dof == 1
    T1.D(1, 1) = 1e-2;
else
    VLinks = T1.VLinks;
    VLinks.Eta = 0.8*VLinks.Eta;
    T1.VLinks = VLinks;

    for i = 1:length(T1.CVRods)
        for j = 1:length(T1.CVRods{i})
            T1.CVRods{1}(1).UpdateAll();
            T1.CVRods{1}(2).UpdateAll();
        end
    end

    % Update Linkage
    T1 = T1.Update();
end


%% Simulation Setup
fs = 1e+3;
tf = 50;
t0 = 0;
t = 0:(1/fs):tf;
N_time = length(t);
n = T1.ndof;

% Repeatable rng
seed = 4;
rng(seed);
q0 = 1.0e-4*randn(n, 1);
qdot0 = zeros(n, 1);
x0 = [q0; qdot0];

% Collocation Object
cf = Collocated_Form(T1);

%% Feasible Target
regen_equilibria = false;
equilibria_dir = fullfile("equilibria", "rcc");
equilibria_file = fullfile(equilibria_dir, "equilibria" + mat_ext);

% Regen Equilibria
if(regen_equilibria || ~exist(equilibria_file, 'file'))
    % Load EquilibriaGVS repo
    addpath(fullfile("..", "GVS-OptimalControl", "EquilibriaGVS"))
    addpath(fullfile("..", "GVS-OptimalControl", "EquilibriaGVS", "functions"))
    
    % Call equilibriaGVS function
    u = zeros(T1.nact, 1);
    equilibria = equilibriaGVS(T1, "input", u);
    % Filter Equilibria
    equilibria = filterEquilibria(equilibria);

    % Save Equilibria
    if(~exist(equilibria_dir, 'dir'))
        mkdir(equilibria_dir);
    end

    save(equilibria_file, "equilibria");
else
    load(equilibria_file);
end

% Stable Equilibrium
q_des = equilibria(:, 1);
q_dot_des = zeros(cf.n, 1);

%% Plot Lyapunov
% Evaluate 
% V = lyapunov(cf, q_des, q_dot_des, q_des, q_dot_des);


function V = lyapunov(cf_obj, q, q_dot, q_des, q_dot_des)
    % Potential Energy at the equilbrium
    [~, Uel_des, Ug_des] = cf_obj.mechanicalEnergy(q_des, q_dot_des);

    % Mechanical Energy at q, qdot
    [T, Uel, Ug] = cf_obj.mechanicalEnergy(q, q_dot);

    % Correction term
    [~, G_theta, K_theta, ~] = cf_obj.transformSystem(q_des, q_dot_des);
    [theta, ~] = cf_obj.transform(q, q_dot);

    % Desired Equilibrium in theta
    [theta_des, ~] = cf_obj.transform(q_des, q_dot_des);

    % Error in collocated variable
    theta_tilde = theta_des - theta;
    theta_tilde_a = theta_tilde(1:cf_obj.m);

    % Compute Lyapunov Function
    V = T + Uel + Ug - Uel_des - Ug_des;
    V = V + (theta_tilde_a')*(G_theta(1:cf_obj.m) + K_theta(1:cf_obj.m));
    % V = V + 0.5*theta_tilde*Kp*theta_tilde;
end