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
q0 = 1.0e-1*randn(n, 1);
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

%% Simulate
% Select Control Law
control_law = "nonlinear_noncollocated_PD_FF";

ODEFUN = @(t, xk) closed_loop(cf, t, xk, "q_des", q_des, ...
                                "control_law", control_law, ...
                                "Kpa", Kpa, "Kda", Kda);

[t_sim, x_sim] = ode15s(ODEFUN, t, x0);
% Column Notation
x_sim = x_sim';

%% Visualization
N_sim = length(t_sim);
z_sim = zeros(2*cf.n, N_sim);

% Convert
for i = 1:N_sim
    z_sim(:, i) = cf.groupTransform(x_sim(:, i));
end

% marker = markers(j);
line_width = 2.5;
line_style = "-";

% Show every simulation
figure
subplot(4, 1, 1)
hold on
plot(t_sim, z_sim(1:cf.m, :), 'LineWidth', line_width, 'LineStyle', line_style)
hold off
grid on
xlabel("$t$ [s]", 'Interpreter', 'latex')
ylabel("$\theta_a$", 'Interpreter', 'latex')
set(gca, 'FontSize', 14)
set(gca, 'GridLineWidth', 1.5)

subplot(4, 1, 2)
hold on
plot(t_sim, z_sim(cf.m + 1:cf.n, :), 'LineWidth', line_width, 'LineStyle', line_style)
hold off
grid on
xlabel("$t$ [s]", 'Interpreter', 'latex')
ylabel("$\theta_u$", 'Interpreter', 'latex')
set(gca, 'FontSize', 14)
set(gca, 'GridLineWidth', 1.5)

subplot(4, 1, 3)
hold on
plot(t_sim, z_sim(cf.n + 1: cf.n + cf.m, :), 'LineWidth', line_width, 'LineStyle', line_style)
hold off
grid on
xlabel("$t$ [s]", 'Interpreter', 'latex')
ylabel("$\dot{\theta}_a$", 'Interpreter', 'latex')
set(gca, 'FontSize', 14)
set(gca, 'GridLineWidth', 1.5)

subplot(4, 1, 4)
hold on
plot(t_sim, z_sim(cf.n + cf.m + 1:end, :), 'LineWidth', line_width, 'LineStyle', line_style)
hold off
grid on
xlabel("$t$ [s]", 'Interpreter', 'latex')
ylabel("$\dot{\theta}_u$", 'Interpreter', 'latex')
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
                            "Kpa", options.Kpa, "Kda", options.Kda);
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
    end

    % Collocated Term
    [u_c, theta_des, theta, theta_dot] = collocated_PD_FF(cf_obj, q_des, q, q_dot, "Kpa", options.Kpa, "Kda", options.Kda);

    % Add noncollocated term
    theta_tilde_u = theta_des(cf_obj.m + 1:end) - theta(cf_obj.m + 1:end);

    %% Compute Nonlinear Gains
    % Potential Hessian (Gamma)
    Gamma = potential_hessian(cf_obj.robot_linkage, q, zeros(cf_obj.n, 1), zeros(cf_obj.m, 1));
    Jh = cf_obj.jacobian(q);
    Gamma_theta = (inv(Jh)')*Gamma*inv(Jh);

    % Compute Damping in collocated variables
    [~, ~, ~, D_theta] = cf_obj.transformSystem(q, zeros(cf_obj.n, 1));

    % Global Gain
    K = 6;

    % Kpu = - 2*Gamma_{a, u} to maximize the convexity
    Kpu = 2.0*K*Gamma_theta(1:cf_obj.m, (cf_obj.m+1):end);
    % Kdu = -2*D_{a, u} to maximize stability
    Kdu = 2.0*K*D_theta(1:cf_obj.m, (cf_obj.m+1):end);

    % Compute Control Law
    u_nc = Kpu*theta_tilde_u - Kdu*theta_dot(cf_obj.m + 1:end);

    % Compose the Actions
    u = u_c + u_nc;
end
