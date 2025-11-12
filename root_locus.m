%%% Hand-made Root Locus %%%
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

%% Linearized System
addpath(fullfile("..", "GVS-OptimalControl", "EquilibriaGVS", "functions"))
[A_lin, B_lin] = linearized_system(T1, q_des, q_dot_des, zeros(T1.nact, 1));

%% Change Coordinates of the linearized system
Jh = cf.jacobian(q_des);
T = blkdiag(Jh, Jh);
% Change of Basis in the Linearized System
A_theta = inv(T)*A_lin*T;
B_theta = inv(T)*B_lin;

%% Show Open-Loop EigenValues
lambda_ol = eig(A_theta);
marker_size = 12;
line_width = 2.0;
% Plot Open-Loop poles on a standard s-plane
% Note: Removed log scale as it's not suitable for negative real parts (stability)
fig_ol = figure;
plot(real(lambda_ol), imag(lambda_ol), 'x', 'MarkerSize', marker_size, 'LineWidth', line_width, "Color", "#de425b")
hold on
grid on
xlabel("Real Part", 'Interpreter', 'latex')
ylabel("Imaginary Part", 'Interpreter', 'latex')
title("Open-Loop Eigenvalues (Poles)", 'Interpreter', 'latex')
set(gca, 'FontSize', 14)
set(gca, 'GridLineWidth', 1.5)
% Add axes at origin
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
hold off

%% Setup for Gain Analysis
% Define the range of gains to test (feel free to change this)
gain_values = 1.0:0.01:10;
num_gains = length(gain_values);

% Get a colormap for plotting
colors = parula(num_gains);

% Define baseline gains (from your original script)
Kpa_base = 1.0*eye(cf.m);
Kda_base = 1.0*eye(cf.m);
Kpu_base = 0.0*ones(cf.m, cf.p);
Kdu_base = 0.0*ones(cf.m, cf.p);

% Plotting parameters
plot_marker_size = 10;
plot_line_width = 2.0;
ol_color = "#de425b"; % Open-loop color from your script

%% 1. Varying Kpa
figure;
hold on
grid on
% Plot Open-Loop poles as reference
plot(real(lambda_ol), imag(lambda_ol), 'x', 'MarkerSize', plot_marker_size, 'LineWidth', plot_line_width, "Color", ol_color);

for i = 1:num_gains
    k_val = gain_values(i);
    
    % Set gains: Vary Kpa, hold others at baseline
    Kpa = k_val * eye(cf.m);
    Kda = Kda_base;
    Kpu = Kpu_base;
    Kdu = Kdu_base;
    
    % Full-State Feedback Gain
    K = [Kpa, Kpu, Kda, Kdu];
    
    % Closed-Loop Eigenvalues
    lambda_cl = eig(A_theta - B_theta*K);
    
    % Plot closed-loop poles
    plot(real(lambda_cl), imag(lambda_cl), 'x', 'MarkerSize', plot_marker_size-2, 'LineWidth', plot_line_width, 'Color', colors(i,:));
end
title('Root Locus: Varying $K_{pa}$', 'Interpreter', 'latex')
xlabel("Real Part", 'Interpreter', 'latex')
ylabel("Imaginary Part", 'Interpreter', 'latex')
set(gca, 'FontSize', 14)
set(gca, 'GridLineWidth', 1.5)
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
% legend(legend_entries_kpa, 'Location', 'best') % Removed as requested
hold off
%% 2. Varying Kda
figure;
hold on
grid on
% Plot Open-Loop poles as reference
plot(real(lambda_ol), imag(lambda_ol), 'x', 'MarkerSize', plot_marker_size, 'LineWidth', plot_line_width, "Color", ol_color);

for i = 1:num_gains
    k_val = gain_values(i);
    
    % Set gains: Vary Kda, hold others at baseline
    Kpa = Kpa_base;
    Kda = k_val * eye(cf.m);
    Kpu = Kpu_base;
    Kdu = Kdu_base;
    
    % Full-State Feedback Gain
    K = [Kpa, Kpu, Kda, Kdu];
    
    % Closed-Loop Eigenvalues
    lambda_cl = eig(A_theta - B_theta*K);
    
    % Plot closed-loop poles
    plot(real(lambda_cl), imag(lambda_cl), 'x', 'MarkerSize', plot_marker_size-2, 'LineWidth', plot_line_width, 'Color', colors(i,:));
end
title('Root Locus: Varying $K_{da}$', 'Interpreter', 'latex')
xlabel("Real Part", 'Interpreter', 'latex')
ylabel("Imaginary Part", 'Interpreter', 'latex')
set(gca, 'FontSize', 14)
set(gca, 'GridLineWidth', 1.5)
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
% legend(legend_entries_kda, 'Location', 'best') % Removed as requested
hold off
%% 3. Varying Kpu
figure;
hold on
grid on
% Plot Open-Loop poles as reference
plot(real(lambda_ol), imag(lambda_ol), 'x', 'MarkerSize', plot_marker_size, 'LineWidth', plot_line_width, "Color", ol_color);

for i = 1:num_gains
    k_val = gain_values(i);
    
    % Set gains: Vary Kpu, hold others at baseline
    Kpa = Kpa_base;
    Kda = Kda_base;
    Kpu = k_val * ones(cf.m, cf.p);
    Kdu = Kdu_base;
    
    % Full-State Feedback Gain
    K = [Kpa, Kpu, Kda, Kdu];
    
    % Closed-Loop Eigenvalues
    lambda_cl = eig(A_theta - B_theta*K);
    
    % Plot closed-loop poles
    plot(real(lambda_cl), imag(lambda_cl), 'x', 'MarkerSize', plot_marker_size-2, 'LineWidth', plot_line_width, 'Color', colors(i,:));
end
title('Root Locus: Varying $K_{pu}$', 'Interpreter', 'latex')
xlabel("Real Part", 'Interpreter', 'latex')
ylabel("Imaginary Part", 'Interpreter', 'latex')
set(gca, 'FontSize', 14)
set(gca, 'GridLineWidth', 1.5)
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
% legend(legend_entries_kpu, 'Location', 'best') % Removed as requested
hold off
%% 4. Varying Kdu
figure;
hold on
grid on
% Plot Open-Loop poles as reference
plot(real(lambda_ol), imag(lambda_ol), 'x', 'MarkerSize', plot_marker_size, 'LineWidth', plot_line_width, "Color", ol_color);

for i = 1:num_gains
    k_val = gain_values(i);
    
    % Set gains: Vary Kdu, hold others at baseline
    Kpa = Kpa_base;
    Kda = Kda_base;
    Kpu = Kpu_base;
    Kdu = k_val * ones(cf.m, cf.p);
    
    % Full-State Feedback Gain
    K = [Kpa, Kpu, Kda, Kdu];
    
    % Closed-Loop Eigenvalues
    lambda_cl = eig(A_theta - B_theta*K);
    
    % Plot closed-loop poles
    plot(real(lambda_cl), imag(lambda_cl), 'x', 'MarkerSize', plot_marker_size-2, 'LineWidth', plot_line_width, 'Color', colors(i,:));
end
title('Root Locus: Varying $K_{du}$', 'Interpreter', 'latex')
xlabel("Real Part", 'Interpreter', 'latex')
ylabel("Imaginary Part", 'Interpreter', 'latex')
set(gca, 'FontSize', 14)
set(gca, 'GridLineWidth', 1.5)
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
% legend(legend_entries_kdu, 'Location', 'best') % Removed as requested
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% END OF MODIFIED SECTION                                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%