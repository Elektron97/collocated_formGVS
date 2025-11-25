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
robot_name = "conical_hsupport";
file_name = "robot_linkage";
mat_ext = ".mat";
% Load Robot and Data
load(fullfile("robots", robot_name, "robot_linkage" + mat_ext));

%% Add DoFs & External Forces
% % Modify DoFs and Modes
% CVRods = T1.CVRods;
% CVRods{1}(2).Phi_dof = [0, 0, 1, 1, 1, 0]';
% CVRods{1}(2).Phi_odr = 2.*CVRods{1}(2).Phi_dof;
% CVRods{1}(2).xi_starfn = @(X) [0, cos(pi/4), sin(pi/4) 1, 0, 0]';
% T1.CVRods = CVRods;

% External Forces
addpath("functions")
T1.CEF = true;

% Update
T1 = T1.Update();

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
% Stable Equilibrium
q_des = equilibria(:, 1);
% Zero velocity
q_dot_des = zeros(cf.n, 1);

%% Linearized System
addpath(fullfile("..", "GVS-OptimalControl", "EquilibriaGVS", "functions"))
[A_lin, B_lin] = linearized_system(T1, q_des, q_dot_des, zeros(T1.nact, 1));

%% Change Coordinates of the linearized system
% Jh = cf.jacobian(q_des);
% T = blkdiag(Jh, Jh);
% % Change of Basis in the Linearized System
% A_theta = inv(T)*A_lin*T;
% B_theta = inv(T)*B_lin;

% Debug
A_theta = A_lin;
B_theta = B_lin;

%% Show Open-Loop EigenValues
lambda_ol = eig(A_theta);
marker_size = 22;
font_size = 28;
grid_linewidth = 3.0;
text_size = 20;
fig = figure;
% Open Loop Plot
plot(real(lambda_ol), imag(lambda_ol), 'x', 'MarkerSize', marker_size, 'LineWidth', 3.0, "Color", red_ol)
hold on
% --- Display Damping Ratio and Natural Frequency (Open Loop) ---
for i = 1:length(lambda_ol)
    val = lambda_ol(i);
    wn = abs(val);
    zeta = -real(val) / (wn + 1e-9); 
    
    txt = sprintf('\\zeta: %.2f\n\\omega_n: %.2f', zeta, wn);
    
    % Offset text positive Y
    text(real(val) + 0.05*abs(real(val)), imag(val) + 0.05*abs(imag(val)), ...
         txt, 'Interpreter', 'tex', 'FontSize', text_size, 'Color', red_ol);
end
% --------------------------------------------------------------
hold off
grid on
xlabel("Real ($\log$ scale)", 'Interpreter', 'latex')
ylabel("Im", 'Interpreter', 'latex')
set(gca, 'XScale', 'log')
set(gca, 'FontSize', font_size)
set(gca, 'GridLineWidth', grid_linewidth)

%% Check on Controllability
for i = 1:length(lambda_ol)
    % Compute HBP test for the open-loop eigenvalues
    pbh = [A_theta - lambda_ol(i)*eye(2*cf.n), B_theta];
    
    % Check the rank
    if(rank(pbh) == 2*cf.n)
        continue;
    else
        disp("The eigenvalue " + num2str(lambda_ol(i)) + " is not controllable.");
    end
end