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
% equilibria_dir = fullfile("equilibria", "rsip");
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

%% Setup for Gain Analysis
% Define baseline gains (from your original script)
Kpa_base = 1.0*eye(cf.m);
Kda_base = 1.0*eye(cf.m);
Kpu_base = 0.0*ones(cf.m, cf.p);
Kdu_base = 0.0*ones(cf.m, cf.p);
% Plotting parameters
plot_marker_size = 16;
plot_line_width = 3.0;
font_size = 20;
grid_linewidth = 2.5;
ol_color = "#de425b"; % Open-loop color from your script

%% 1. Varying Kpa
% Define the range of gains for Kpa
gain_values_kpa = 1.0:0.1:10; % <-- CHANGE THIS RANGE
num_gains = length(gain_values_kpa);
% colors = parula(num_gains);
colors = SoFFTColormap(num_gains);
figure;
hold on
grid on
% Plot Open-Loop poles as reference
plot(real(lambda_ol), imag(lambda_ol), 'x', 'MarkerSize', plot_marker_size, 'LineWidth', plot_line_width, "Color", ol_color);
for i = 1:num_gains
    k_val = gain_values_kpa(i);
    
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
set(gca, 'XScale', 'log')
set(gca, 'FontSize', font_size)
set(gca, 'GridLineWidth', grid_linewidth)
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

% --- COLORBAR UPDATE ---
colormap(colors); % Set the colormap to match the colors used
c = colorbar;
clim([min(gain_values_kpa) max(gain_values_kpa)]); % Scale colorbar to gain range
c.Label.String = '$K_{pa}$ Magnitude';
c.Label.Interpreter = 'latex';
c.Label.FontSize = font_size;
hold off

%% 2. Varying Kda
% Define the range of gains for Kda
gain_values_kda = 1.0:0.1:10.0; % <-- CHANGE THIS RANGE
num_gains = length(gain_values_kda);
colors = SoFFTColormap(num_gains);
figure;
hold on
grid on
% Plot Open-Loop poles as reference
plot(real(lambda_ol), imag(lambda_ol), 'x', 'MarkerSize', plot_marker_size, 'LineWidth', plot_line_width, "Color", ol_color);
for i = 1:num_gains
    k_val = gain_values_kda(i);
    
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
set(gca, 'FontSize', font_size)
set(gca, 'GridLineWidth', grid_linewidth)
set(gca, 'XScale', 'log')
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

% --- COLORBAR UPDATE ---
colormap(colors);
c = colorbar;
clim([min(gain_values_kda) max(gain_values_kda)]); % Scale colorbar
c.Label.String = '$K_{da}$ Magnitude';
c.Label.Interpreter = 'latex';
c.Label.FontSize = font_size;
hold off

%% 3. Varying Kpu
% Define the range of gains for Kpu
gain_values_kpu = -5.0:0.1:0.0; % <-- CHANGE THIS RANGE
num_gains = length(gain_values_kpu);
colors = SoFFTColormap(num_gains);
figure;
hold on
grid on
% Plot Open-Loop poles as reference
plot(real(lambda_ol), imag(lambda_ol), 'x', 'MarkerSize', plot_marker_size, 'LineWidth', plot_line_width, "Color", ol_color);
for i = 1:num_gains
    k_val = gain_values_kpu(i);
    
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
set(gca, 'FontSize', font_size)
set(gca, 'GridLineWidth', grid_linewidth)
set(gca, 'XScale', 'log')
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

% --- COLORBAR UPDATE ---
colormap(colors);
c = colorbar;
clim([min(gain_values_kpu) max(gain_values_kpu)]); % Scale colorbar
c.Label.String = '$K_{pu}$ Magnitude';
c.Label.Interpreter = 'latex';
c.Label.FontSize = font_size;
hold off

%% 4. Varying Kdu
% Define the range of gains for Kdu
gain_values_kdu = 0.0:0.01:0.35; % <-- CHANGE THIS RANGE
num_gains = length(gain_values_kdu);
colors = SoFFTColormap(num_gains);
figure;
hold on
grid on
% Plot Open-Loop poles as reference
plot(real(lambda_ol), imag(lambda_ol), 'x', 'MarkerSize', plot_marker_size, 'LineWidth', plot_line_width, "Color", ol_color);
for i = 1:num_gains
    k_val = gain_values_kdu(i);
    
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
set(gca, 'FontSize', font_size)
set(gca, 'GridLineWidth', grid_linewidth)
set(gca, 'XScale', 'log')
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

% --- COLORBAR UPDATE ---
colormap(colors);
c = colorbar;
clim([min(gain_values_kdu) max(gain_values_kdu)]); % Scale colorbar
c.Label.String = '$K_{du}$ Magnitude';
c.Label.Interpreter = 'latex';
c.Label.FontSize = font_size;
hold off

%% ColorMap
function map = SoFFTColormap(n)
    % 1. Load the fixed 256x3 data
    baseMap = load("SoFFTColormap.mat");    
    m = size(baseMap.SoFFTColormap, 1);
    
    % 2. Handle input arguments (Mimic parula's behavior)
    if nargin < 1
        f = get(groot, 'CurrentFigure');
        if isempty(f)
            n = m; % Default to 256 if no figure exists
        else
            n = size(f.Colormap, 1); % Default to current figure's colormap length
        end
    end
    
    % 3. Interpolate from 256 down (or up) to n
    map = interp1(1:m, baseMap.SoFFTColormap, linspace(1, m, n), 'linear');
end