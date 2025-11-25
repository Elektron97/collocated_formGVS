%%% Collocated Form Control %%%
clear all
close all
clc

%% Load and Startup SoRoSim
% Clean StartUp
diff_sorosim_path = fullfile("SoRoSim", "Differentiable_SoRoSim");
if exist(diff_sorosim_path, 'dir')
    cd(diff_sorosim_path)
    startup
    % Switch again to the current directory
    [current_path, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
    cd(current_path)
end

%% Load Data
% robot_name = "rsip";
% robot_name = "rsip_extreme";
% robot_name = "conical_hsupport";
robot_name = "cable";
file_name = "robot_linkage";
mat_ext = ".mat";
% Load Robot and Data
load(fullfile("robots", robot_name, "robot_linkage" + mat_ext));

%% Add DoFs & External Forces
% % Modify DoFs and Modes
% CVRods = T1.CVRods;
% CVRods{1}(2).Phi_dof = [0, 0, 1, 1, 1, 0]';
% CVRods{1}(2).Phi_odr = 2.*CVRods{1}(2).Phi_dof;
% T1.CVRods = CVRods;
% 
% % External Forces
% addpath("functions")
% T1.CEF = true;
% T1 = T1.Update();

% % Modify Young Modulus
% VLinks = T1.VLinks;
% VLinks.E = 1.0e-1*VLinks.E;
% T1.VLinks = VLinks;
% T1.CVRods{1}(2).UpdateMEG();
% T1 = T1.Update();

% Joint Damping
T1.D(1, 1) = 0.5;

%% Colors
blue_sofft   = "#086788";   % Non-Collocated (Full)
red_target   = "#f06543";
grey_mid     = "#353B45";
red_ol       = "#de425b";   % Open Loop (kept for Eigenvalues plot only)
green_col    = "#00a066";   % Collocated Only
yellow_nokpu = "#ffa600";   % Non-Collocated (Kpu = 0)

% Create Figure Path
figure_path = fullfile("figures", robot_name);
% Create if not exists
if(~exist(figure_path, 'dir'))
    mkdir(figure_path);
end

%% Simulation Setup
fs = 1e+3;
tf = 10;
t0 = 0;
t = 0:(1/fs):tf;
N_time = length(t);
n = T1.ndof;

% Repeatable rng
seed = 4;
rng(seed);
q0 = 0.0e+0*randn(n, 1);
qdot0 = zeros(n, 1);
x0 = [q0; qdot0];

% Collocation Object
cf = Collocated_Form(T1);

%% Feasible Target
% Select Equilibria for input    
u = 0.0.*ones(T1.nact, 1);
regen_equilibria = true;
equilibria_dir = fullfile("equilibria", robot_name);
equilibria_file = fullfile(equilibria_dir, "equilibria" + "_" + num2str(u) + mat_ext);
% Regen Equilibria
if(regen_equilibria || ~exist(equilibria_file, 'file'))
    % Load EquilibriaGVS repo
    addpath(fullfile("..", "GVS-OptimalControl", "EquilibriaGVS"))
    addpath(fullfile("..", "GVS-OptimalControl", "EquilibriaGVS", "functions"))
    
    % Call equilibriaGVS function
    equilibria = equilibriaGVS(T1, "input", u);
    % Filter Equilibria
    equilibria = filterEquilibria(equilibria, "angle_mask", [true, false, false]);
    % Save Equilibria
    if(~exist(equilibria_dir, 'dir'))
        mkdir(equilibria_dir);
    end
    % save(equilibria_file, "equilibria");
else
    load(equilibria_file);
end

% true: stable equilibrium | false: unstable equilibrium
is_stable = false;

if is_stable
    % Stable Equilibrium
    q_des = equilibria(:, 1);
else
    % Unstable Equilibrium
    q_des = equilibria(:, 2);
end
% Zero velocity
q_dot_des = zeros(cf.n, 1);

%% Linearized System
addpath(fullfile("..", "GVS-OptimalControl", "EquilibriaGVS", "functions"))
[A_lin, B_lin] = linearized_system(T1, q_des, q_dot_des, zeros(T1.nact, 1));

%% Change Coordinates of the linearized system
Jh = cf.jacobian(q_des);
T = blkdiag(Jh, Jh);
% Change of Basis in the Linearized System
A_theta = inv(T)*A_lin*T;
B_theta = inv(T)*B_lin;