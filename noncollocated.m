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
q0 = 0.1*randn(n, 1);
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

%% Linearized System
addpath(fullfile("..", "GVS-OptimalControl", "EquilibriaGVS", "functions"))
[A_lin, B_lin] = linearized_system(T1, q_des, q_dot_des, zeros(T1.nact, 1));

%% Design LQR
% Q = blkdiag(eye(cf.m), eye(cf.p), eye(cf.m), eye(cf.p));
% R = 1e+0*eye(cf.m);
% 
% % Solve Riccati
% [K_lqr, ~, lambda_cl] = lqr(A_lin, B_lin, Q, R);
% 
% % Extract Gains
% Kpa = K_lqr(1:cf.m);
% Kpu = K_lqr(cf.m + 1:cf.n);
% Kda = K_lqr(cf.n + 1: cf.n + cf.m);
% Kdu = K_lqr(cf.n + cf.m + 1:end);

% Numerical Gains
Kpa = 1e+2;
Kda = 1e+1;
Kpu = 2e+0*[1.0, 1.0];
Kdu = 5e+0*[1.0, 0.0];

% Trials
Kpas = {Kpa*eye(cf.m), Kpa*eye(cf.m), Kpa*eye(cf.m), Kpa*eye(cf.m)};
Kdas = {Kda*eye(cf.m), Kda*eye(cf.m), Kda*eye(cf.m), Kda*eye(cf.m)};
% Kpus = {zeros(cf.m, cf.p), zeros(cf.m, cf.p),   zeros(cf.m, cf.p),  zeros(cf.m, cf.p)};
% Kdus = {zeros(cf.m, cf.p), zeros(cf.m, cf.p),   zeros(cf.m, cf.p),  zeros(cf.m, cf.p)};
% Kpus = {zeros(cf.m, cf.p), zeros(cf.m, cf.p),   zeros(cf.m, cf.p),  zeros(cf.m, cf.p)};
% Kdus = {zeros(cf.m, cf.p), 0.5*Kdu,             Kdu,                1.2*Kdu};
% Kpus = {zeros(cf.m, cf.p), 0.5*Kpu,             Kpu,                1.2*Kpu};
% Kdus = {zeros(cf.m, cf.p), zeros(cf.m, cf.p),   zeros(cf.m, cf.p),  zeros(cf.m, cf.p)};
Kpus = {zeros(cf.m, cf.p), 0.5*Kpu,             Kpu,                1.2*Kpu};
Kdus = {zeros(cf.m, cf.p), 0.5*Kdu,             Kdu,                1.2*Kdu};


%% Compute Closed-Loop Eigenvalues (Numerical)
% Open Loop EigenValues
lambda_ol = eig(A_lin);

% Closed-Loop EigenValues
[M, ~, ~, ~] = cf.dynamicMatrices(q_des, q_dot_des);
Jh = cf.jacobian(q_des);

% Analyze all the cases
for i = 1:length(Kdus)
    Kp = (Jh')*[Kpas{i}, Kpus{i}; zeros(cf.p, cf.n)]*Jh;
    Kd = (Jh')*[Kdas{i}, Kdus{i}; zeros(cf.p, cf.n)]*Jh;
    
    A_cl = A_lin + [zeros(cf.n, cf.n), zeros(cf.n, cf.n); -M\Kp, -M\Kd];
    lambda_cl{i} = eig(A_cl);
end

% Palette
poles_palette = [blue_sofft, "#00878d", "#00a065", "#8cae22", "#ffa600"];

figure
% Open Loop
plot(real(lambda_ol), imag(lambda_ol), 'x', 'MarkerSize', 10, 'LineWidth', 2.0, "Color", blue_sofft)
hold on
for i = 1:length(lambda_ol)
    label = num2str(i);
    text(real(lambda_ol(i)) + 0.2, imag(lambda_ol(i)) + 0.2, label, 'FontSize', 16, 'Color', blue_sofft, 'FontWeight', 'bold');
end

% Closed Loop
for i = 1:length(Kdus)
    plot(real(lambda_cl{i}), imag(lambda_cl{i}), 'x', 'MarkerSize', 10, 'LineWidth', 2.0, "Color", poles_palette(i + 1))
end
hold off
grid on
xlabel("Real")
ylabel("Im")
set(gca, 'XScale', 'log')

%% Simulate
% Store Controllers results
t_sim = cell(length(Kpus), 1);
x_sim = cell(length(Kpus), 1);

% Select Control Law
control_law = "noncollocated_PD_FF";

% for i = 1
for i = 1:length(Kpus)
    disp("Simulation n.: " + num2str(i))
    ODEFUN = @(t, xk) closed_loop(cf, t, xk, "q_des", q_des, "control_law", control_law, ...
                                        "Kpa", Kpas{i}, "Kda", Kdas{i}, ...
                                        "Kpu", Kpus{i}, "Kdu", Kdus{i});

    [t_sim{i}, x_sim{i}] = ode15s(ODEFUN, t, x0);
    % Column Notation
    x_sim{i} = x_sim{i}';
end

%% Visualization
% T1.plotqt(t_sim{1}, x_sim{1}', "record", false)

%% Convert in Collocated Variables
z_sim = cell(length(Kpus), 1);

% Target
[theta_des, theta_dot_des] = cf.transform(q_des, q_dot_des);

% Line Styles
% markers = ["o", "+", "*", "x", "square"];
% line_styles = ["-", "--", ":", "-."];

% Active and Passive Colors
% active_colors = ["#de425b", "#ed7883", "#f9a7ac", "#ffd5d7"];
% passive_colors = ["#086888", "#5191ae", "#5191ae", "#bae9ff"];
% active_colors = ["#086888", "#009084", "#67ab3b", "#ffa600"];
% passive_colors = active_colors;
% palette = [blue_sofft, "#009084", "#67ab3b", "#ffa600"];
palette = poles_palette(2:end);

% Show
figure

% Convert in theta
for j = 1:length(Kpus)
    N_sim = length(t_sim{j});
    z_sim{j} = zeros(2*cf.n, N_sim);

    % Convert
    for i = 1:N_sim
        z_sim{j}(:, i) = cf.groupTransform(x_sim{j}(:, i));
    end

    % marker = markers(j);
    line_width = 2.5;
    line_style = "-";

    % Show every simulation
    subplot(4, 1, 1)
    hold on
    plot(t_sim{j}, z_sim{j}(1:cf.m, :), 'LineWidth', line_width, 'LineStyle', line_style, 'Color', palette(j))
    hold off
    grid on
    xlabel("$t$ [s]", 'Interpreter', 'latex')
    ylabel("$\theta_a$", 'Interpreter', 'latex')

    subplot(4, 1, 2)
    hold on
    plot(t_sim{j}, z_sim{j}(cf.m + 1:cf.n, :), 'LineWidth', line_width, 'LineStyle', line_style, 'Color', palette(j))
    hold off
    grid on
    xlabel("$t$ [s]", 'Interpreter', 'latex')
    ylabel("$\theta_u$", 'Interpreter', 'latex')

    subplot(4, 1, 3)
    hold on
    plot(t_sim{j}, z_sim{j}(cf.n + 1: cf.n + cf.m, :), 'LineWidth', line_width, 'LineStyle', line_style, 'Color', palette(j))
    hold off
    grid on
    xlabel("$t$ [s]", 'Interpreter', 'latex')
    ylabel("$\dot{\theta}_a$", 'Interpreter', 'latex')

    subplot(4, 1, 4)
    hold on
    plot(t_sim{j}, z_sim{j}(cf.n + cf.m + 1:end, :), 'LineWidth', line_width, 'LineStyle', line_style, 'Color', palette(j))
    hold off
    grid on
    xlabel("$t$ [s]", 'Interpreter', 'latex')
    ylabel("$\dot{\theta}_u$", 'Interpreter', 'latex')
end

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

%% Build Rectangular Gains
function A = build_rectangular_gains(singular_values, m, n)
    %generateMatrixFromSVs Generates an m-by-n matrix with specified singular values.

    %% 1. Input Validation
    if any(singular_values < 0)
        error('Singular values must be non-negative.');
    end
    
    k = length(singular_values);
    if k > min(m, n)
        error('The number of singular values cannot exceed min(m, n).');
    end

    %% 2. Construct the Sigma Matrix
    Sigma = zeros(m, n);
    Sigma(1:k, 1:k) = diag(singular_values);

    %% 3. Generate Random Orthogonal Matrices U and V
    % The Q factor from the QR decomposition of a random matrix is a 
    % random orthogonal matrix.
    [U, ~] = qr(randn(m, m), 0); % Using the "economy" size QR
    [V, ~] = qr(randn(n, n), 0); % Using the "economy" size QR

    %% 4. Assemble the Final Matrix
    A = U * Sigma * V';

end
