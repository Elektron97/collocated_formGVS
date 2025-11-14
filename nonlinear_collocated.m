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
addpath("functions")

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
q0 = 1.0e+0*randn(n, 1);
qdot0 = zeros(n, 1);
x0 = [q0; qdot0];

% Collocation Object
cf = Collocated_Form(T1);

%% Feasible Target
regen_equilibria = false;
equilibria_dir = fullfile("equilibria", robot_name);
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

%% Gains
Kpa = eye(cf.m);
Kda = eye(cf.m);

%% Simulate for multiple K values
% Define K values to test
% K_values = linspace(0, 50, 8);
K_values = 0:0.5:50;
N_sims = length(K_values);

% Get control law
control_law = "nonlinear_noncollocated_PD_FF";

% Pre-allocate cell arrays to store results
t_sim_all = cell(N_sims, 1);
x_sim_all = cell(N_sims, 1);
z_sim_all = cell(N_sims, 1);

% Pre-calculate desired collocated position vector
[theta_des, theta_dot_des] = cf.transform(q_des, q_dot_des);
% Pre-allocate array for ISE values
ISE_values = zeros(N_sims, 1);

fprintf('Running %d simulations for different K values...\n', N_sims);

% Loop over each K value
for k = 1:N_sims
    K_val = K_values(k);
    K_mat = K_val * eye(cf.p);
    
    fprintf('  Simulating with K = %.2f\n', K_val);
    
    % Define the ODE function for the current K
    ODEFUN = @(t, xk) closed_loop(cf, t, xk, "q_des", q_des, ...
                                    "control_law", control_law, ...
                                    "Kpa", Kpa, "Kda", Kda, ...
                                    "K", K_mat); % Use the current K
                                
    % Run the simulation
    [t_sim, x_sim] = ode15s(ODEFUN, t, x0);
    
    % Store results (transpose x_sim for consistency with original)
    t_sim_all{k} = t_sim;
    x_sim_all{k} = x_sim';
    
    % Post-process: Convert to z_sim
    N_sim_steps = length(t_sim);
    x_sim_k_transposed = x_sim_all{k}; % Already transposed
    z_sim_k = zeros(2*cf.n, N_sim_steps);
    for i = 1:N_sim_steps
        z_sim_k(:, i) = cf.groupTransform(x_sim_k_transposed(:, i));
    end
    z_sim_all{k} = z_sim_k;
    
    % Get time and state for this simulation
    t_sim_k = t_sim_all{k};         % [N_sim_steps x 1]
    z_sim_k = z_sim_all{k};         % [2*n x N_sim_steps]
    
    % Error
    error_k = z_sim_k - [theta_des; theta_dot_des];
    % error_k = z_sim_k(1:cf.n, :) - theta_des;
    % error_k = z_sim_k(cf.n + 1:end, :) - theta_dot_des;
    
    % Calculate the squared norm of the error at each time step: e(t)'*e(t)
    % sum(error_k.^2, 1) squares elements and sums down columns
    squared_error_norm_k = sum(error_k.^2, 1); % [1 x N_sim_steps]
    
    % Integrate over time using trapezoidal rule
    % Note: t_sim_k must be [N_sim_steps x 1] and squared_error_norm_k must be [1 x N_sim_steps] or [N_sim_steps x 1]
    ISE_values(k) = trapz(t_sim_k, squared_error_norm_k);
end

fprintf('Simulations complete. Plotting results...\n');

%% Visualization
% Load Custom Colormap
load("SoFFTColormap.mat")

% Get colormap
colors = turbo(N_sims);
line_width = 2.0;
line_style = "-";

% Create legend entries
legend_entries = cell(N_sims, 1);
for k = 1:N_sims
    legend_entries{k} = sprintf('K = %.2f', K_values(k));
end

% Create the figure
figure

% Subplot 1: \theta_a
subplot(4, 1, 1)
hold on
for k = 1:N_sims
    plot(t_sim_all{k}, z_sim_all{k}(1:cf.m, :), ...
        'LineWidth', line_width, 'LineStyle', line_style, 'Color', colors(k, :))
end
hold off
grid on
xlabel("$t$ [s]", 'Interpreter', 'latex')
ylabel("$\theta_a$", 'Interpreter', 'latex')
set(gca, 'FontSize', 14)
set(gca, 'GridLineWidth', 1.5)
% Add legend to the first plot
legend(legend_entries, 'Location', 'best')

% Subplot 2: \theta_u
subplot(4, 1, 2)
hold on
for k = 1:N_sims
    plot(t_sim_all{k}, z_sim_all{k}(cf.m + 1:cf.n, :), ...
        'LineWidth', line_width, 'LineStyle', line_style, 'Color', colors(k, :))
end
hold off
grid on
xlabel("$t$ [s]", 'Interpreter', 'latex')
ylabel("$\theta_u$", 'Interpreter', 'latex')
set(gca, 'FontSize', 14)
set(gca, 'GridLineWidth', 1.5)

% Subplot 3: \dot{\theta}_a
subplot(4, 1, 3)
hold on
for k = 1:N_sims
    plot(t_sim_all{k}, z_sim_all{k}(cf.n + 1: cf.n + cf.m, :), ...
        'LineWidth', line_width, 'LineStyle', line_style, 'Color', colors(k, :))
end
hold off
grid on
xlabel("$t$ [s]", 'Interpreter', 'latex')
ylabel("$\dot{\theta}_a$", 'Interpreter', 'latex')
set(gca, 'FontSize', 14)
set(gca, 'GridLineWidth', 1.5)

% Subplot 4: \dot{\theta}_u
subplot(4, 1, 4)
hold on
for k = 1:N_sims
    plot(t_sim_all{k}, z_sim_all{k}(cf.n + cf.m + 1:end, :), ...
        'LineWidth', line_width, 'LineStyle', line_style, 'Color', colors(k, :))
end
hold off
grid on
xlabel("$t$ [s]", 'Interpreter', 'latex')
ylabel("$\dot{\theta}_u$", 'Interpreter', 'latex')
set(gca, 'FontSize', 14)
set(gca, 'GridLineWidth', 1.5)

%% Plot ISE vs. K
figure
plot(K_values, ISE_values, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'Color', red_target)
hold on
% Find the minimum ISE
[min_ISE, min_idx] = min(ISE_values);
min_K = K_values(min_idx);
% Plot a marker for the minimum
plot(min_K, min_ISE, 'p', 'MarkerSize', 14, 'MarkerFaceColor', blue_sofft, 'MarkerEdgeColor', 'k')
legend('ISE', sprintf('Min ISE at K = %.2f', min_K), 'Location', 'best')

xlabel('Gain $K$', 'Interpreter', 'latex')
ylabel('ISE: $\int ||e(t)||^2 dt$', 'Interpreter', 'latex')
title('ISE (Integral of Squared Error) vs. Gain $K$', 'Interpreter', 'latex')
grid on
set(gca, 'FontSize', 14)
set(gca, 'GridLineWidth', 1.5)

%% Functions
function x_dot = dynamics(robot_linkage, t, x, u)
    % dynamicsSolver function
    [y,~,~,~] = robot_linkage.dynamicsSolver(t, x, u);
    
    % dxdt = [q_dot; q_2dot]
    x_dot = [x(robot_linkage.ndof + 1:end); y(1:robot_linkage.ndof)];
end
function x_dot = closed_loop(cf, t, x, options)
    arguments
        cf; 
        t; 
        x;
        options.control_law = "collocated_PD_FF";
        options.Kpa = eye(cf.m);
        options.Kda = eye(cf.m);
        options.Kpu = zeros(cf.m, cf.p);
        options.Kdu = zeros(cf.m, cf.p);
        options.K = eye(cf.p);
        options.q_des = zeros(2*cf.n, 1);
    end
    % Compute Control Action
    u = 0.0*ones(cf.m, 1);
    % Select control law
    switch options.control_law
        case "autonomous"
            u = zeros(cf.m, 1);
        case "collocated_PD_FF"
            [u, ~, ~, ~] = collocated_PD_FF(cf, options.q_des, x(1:cf.n), x(cf.n + 1:end), "Kpa", options.Kpa, "Kda", options.Kda);
        case "noncollocated_PD_FF"
            u = noncollocated_PD_FF(cf, options.q_des, x(1:cf.n), x(cf.n + 1:end), ...
                                        "Kpa", options.Kpa, "Kda", options.Kda, ...
                                        "Kpu", options.Kpu, "Kdu", options.Kdu);
        case "nonlinear_noncollocated_PD_FF"
            u = nonlinear_noncollocated_PD_FF(cf, options.q_des, x(1:cf.n), x(cf.n + 1:end), ...
                            "Kpa", options.Kpa, "Kda", options.Kda, ...
                            "K", options.K);
        otherwise
            warning("Control Law not supported.")
    end
    % Compute Dynamics
    x_dot = dynamics(cf.robot_linkage, t, x, u);
end

%% Collocated Law
function [u, theta_des, theta, theta_dot] = collocated_PD_FF(cf_obj, q_des, q, q_dot, options)
    arguments
        % Collocation Form object: useful for conversion
        cf_obj
        % Desired state (in joint space, not in the collocated variables)
        q_des
        % Feedback terms
        q
        q_dot
        % Gains
        options.Kpa = eye(cf_obj.m);
        options.Kda = eye(cf_obj.m);
    end
    %% Compute Compensation Terms
    [~, G_theta_des, K_theta_des, ~] = cf_obj.transformSystem(q_des, zeros(cf_obj.n, 1));
    % Extract collocated parts
    Ga = G_theta_des(1:cf_obj.m);
    Ka = K_theta_des(1:cf_obj.m);
    %% Convert Desired Configuration in Collocated Form
    [theta_des, ~] = cf_obj.transform(q_des, zeros(cf_obj.n, 1));
    theta_des_a = theta_des(1:cf_obj.m);
    %% Convert Feedback in the Collocated Variables
    [theta, theta_dot] = cf_obj.transform(q, q_dot);
    % Extract Collocated part
    theta_a = theta(1:cf_obj.m);
    theta_dot_a = theta_dot(1:cf_obj.m);
    %% Compute Control Law
    u = options.Kpa*(theta_des_a - theta_a) - options.Kda*(theta_dot_a) + Ga + Ka;
end

%% Non-collocated Law
function u = noncollocated_PD_FF(cf_obj, q_des, q, q_dot, options)
    arguments
        % Collocation Form object: useful for conversion
        cf_obj
        % Desired state (in joint space, not in the collocated variables)
        q_des
        % Feedback terms
        q
        q_dot
        % Gains
        options.Kpa = eye(cf_obj.m);
        options.Kda = eye(cf_obj.m);
        options.Kpu = eye(cf_obj.p);
        options.Kdu = eye(cf_obj.p);
    end
    % Collocated Term
    [u_c, theta_des, theta, theta_dot] = collocated_PD_FF(cf_obj, q_des, q, q_dot, "Kpa", options.Kpa, "Kda", options.Kda);
    % Add noncollocated term
    theta_tilde_u = theta_des(cf_obj.m + 1:end) - theta(cf_obj.m + 1:end);
    u_nc = options.Kpu*theta_tilde_u - options.Kdu*theta_dot(cf_obj.m + 1:end);
    % Compose the Actions
    u = u_c + u_nc;
end
function u = nonlinear_noncollocated_PD_FF(cf_obj, q_des, q, q_dot, options)
    arguments
        % Collocation Form object: useful for conversion
        cf_obj
        % Desired state (in joint space, not in the collocated variables)
        q_des
        % Feedback terms
        q
        q_dot
        % Gains
        options.Kpa = eye(cf_obj.m);
        options.Kda = eye(cf_obj.m);
        options.K = eye(cf_obj.p);
    end

    % Collocated Term
    [u_c, theta_des, theta, theta_dot] = collocated_PD_FF(cf_obj, q_des, q, q_dot, "Kpa", options.Kpa, "Kda", options.Kda);

    % Add noncollocated term
    theta_tilde_u = theta_des(cf_obj.m + 1:end) - theta(cf_obj.m + 1:end);

    %% Compute Nonlinear Gains
    % Potential Hessian (Gamma)
    % Gamma = potential_hessian(cf_obj.robot_linkage, q, zeros(cf_obj.n, 1), zeros(cf_obj.m, 1));

    % Linear (Constant Gamma) case
    Gamma = potential_hessian(cf_obj.robot_linkage, q_des, zeros(cf_obj.n, 1), zeros(cf_obj.m, 1));
    Jh = cf_obj.jacobian(q);
    Gamma_theta = (inv(Jh)')*Gamma*inv(Jh);

    % Compute Damping in collocated variables
    % [~, ~, ~, D_theta] = cf_obj.transformSystem(q, q_dot);
    [~, ~, ~, D_theta] = cf_obj.transformSystem(q_des, zeros(cf_obj.n, 1));

    % Kpu = - 2*Gamma_{a, u} to maximize the convexity
    Kpu = -2.0*Gamma_theta(1:cf_obj.m, (cf_obj.m+1):end)*options.K;
    
    % Kdu = -2*D_{a, u} to maximize stability
    Kdu = -2.0*D_theta(1:cf_obj.m, (cf_obj.m+1):end)*options.K;
    % Compute Control Law
    u_nc = Kpu*theta_tilde_u - Kdu*theta_dot(cf_obj.m + 1:end);
    % Compose the Actions
    u = u_c + u_nc;
end