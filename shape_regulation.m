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
robot_name = "rsip";
% robot_name = "rsip_extreme";
% robot_name = "conical_hsupport";

file_name = "robot_linkage";
mat_ext = ".mat";

% Load Robot and Data
load(fullfile("robots", robot_name, "robot_linkage" + mat_ext));

%% Add DoFs & External Forces
% Modify DoFs and Modes
CVRods = T1.CVRods;
CVRods{1}(2).Phi_dof = [0, 0, 1, 1, 1, 0]';
CVRods{1}(2).Phi_odr = 2.*CVRods{1}(2).Phi_dof;
T1.CVRods = CVRods;

% External Forces
addpath("functions")
T1.CEF = true;

T1 = T1.Update();
T1.D(1, 1) = 0.5;

%% Colors
blue_sofft  = "#086788";   % Non-Collocated (Full)
red_target  = "#f06543";
grey_mid    = "#858583";
red_ol      = "#de425b";     % Open Loop
green_col   = "#2ca02c";    % Collocated Only

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
regen_equilibria = true;
equilibria_dir = fullfile("equilibria", robot_name);
equilibria_file = fullfile(equilibria_dir, "equilibria" + mat_ext);

% Regen Equilibria
if(regen_equilibria || ~exist(equilibria_file, 'file'))
    % Load EquilibriaGVS repo
    addpath(fullfile("..", "GVS-OptimalControl", "EquilibriaGVS"))
    addpath(fullfile("..", "GVS-OptimalControl", "EquilibriaGVS", "functions"))
    
    % Call equilibriaGVS function
    u = 0.0.*ones(T1.nact, 1);
    equilibria = equilibriaGVS(T1, "input", u);
    % Filter Equilibria
    equilibria = filterEquilibria(equilibria, "angle_mask", [true, false, false]);
    % Save Equilibria
    if(~exist(equilibria_dir, 'dir'))
        mkdir(equilibria_dir);
    end
    save(equilibria_file, "equilibria");
else
    load(equilibria_file);
end

% % Stable Equilibrium
% q_des = equilibria(:, 1);
% Unstable Equilibrium
q_des = equilibria(:, 2);
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
marker_size = 20;
font_size = 20;
grid_linewidth = 2.5;
text_size = 18;

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
% set(gca, 'XScale', 'log')
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

%% Pole Placement or LQR
linear_control = "lqr"; % "pole_placement", "lqr"
lambda_cl = lambda_ol;

switch(linear_control)
    case "pole_placement"
        poles = [-1, -2, -3, -4, -5, -6];
        [K_pp, precision] = place(A_theta, B_theta, poles);
        % Extract Gains
        Kpa = K_pp(1:cf.m);
        Kpu = K_pp(cf.m + 1:cf.n);
        Kda = K_pp(cf.n + 1: cf.n + cf.m);
        Kdu = K_pp(cf.n + cf.m + 1:end);
        % Compute Closed-Loop Eigenvalues
        lambda_cl = eig(A_theta - B_theta*K_pp);
        
    case "lqr"
        %% Apply LQR
        Q = blkdiag(1e+3*eye(cf.m), 1e+0*eye(cf.p), 1e+3*eye(cf.m), 1e+0*eye(cf.p));
        R = 1e+0*eye(cf.m);
        
        % Solve Riccati
        [K_lqr, P, lambda_cl] = lqr(A_theta, B_theta, Q, R);
        
        % Extract Gains
        Kpa = K_lqr(1:cf.m);
        Kpu = K_lqr(cf.m + 1:cf.n);
        Kda = K_lqr(cf.n + 1: cf.n + cf.m);
        Kdu = K_lqr(cf.n + cf.m + 1:end);
    otherwise
        error("Linear Controller not supported.")
end

% Plot Closed-Loop Eigenvalues
figure(fig)
hold on
plot(real(lambda_cl), imag(lambda_cl), 'x', 'MarkerSize', marker_size, 'LineWidth', 3.0, "Color", blue_sofft)
% --- Display Damping Ratio and Natural Frequency (Closed Loop) ---
for i = 1:length(lambda_cl)
    val = lambda_cl(i);
    wn = abs(val);
    zeta = -real(val) / (wn + 1e-9); 
    
    txt = sprintf('\\zeta: %.2f\n\\omega_n: %.2f', zeta, wn);
    
    % Offset text negative Y (to distinguish from open loop if close)
    text(real(val) + 0.1*abs(real(val)), imag(val) - 0.1*abs(imag(val)), ...
         txt, 'Interpreter', 'tex', 'FontSize', text_size, 'Color', blue_sofft);
end
hold off
legend("Open-Loop", "Closed-Loop")

%% Simulate
% 1. Simulate Open Loop (Zero Gains)
ODEFUN_OL = @(t, xk) closed_loop(cf, t, xk, "control_law", "autonomous");
[t_sim_ol, x_sim_ol] = ode15s(ODEFUN_OL, t, x0);
x_sim_ol = x_sim_ol';

% 2. Simulate Non-Collocated (Full Gains: Kpa, Kda, Kpu, Kdu)
control_law = "noncollocated_PD_FF";
ODEFUN_FULL = @(t, xk) closed_loop(cf, t, xk, "q_des", q_des, ...
                                "control_law", control_law, ...
                                "Kpa", Kpa, "Kda", Kda, ...
                                "Kpu", Kpu, "Kdu", Kdu);
[t_sim, x_sim] = ode15s(ODEFUN_FULL, t, x0);
x_sim = x_sim';

% 3. Simulate Collocated Only (Partial Gains: Kpa, Kda only)
control_law_col = "collocated_PD_FF";
ODEFUN_COL = @(t, xk) closed_loop(cf, t, xk, "q_des", q_des, ...
                                "control_law", control_law_col, ...
                                "Kpa", Kpa, "Kda", Kda); % Kpu and Kdu are ignored by this law
[t_sim_col, x_sim_col] = ode15s(ODEFUN_COL, t, x0);
x_sim_col = x_sim_col';

%% Visualization
% --- Process Non-Collocated Data ---
N_sim = length(t_sim);
z_sim = zeros(2*cf.n, N_sim);
for i = 1:N_sim
    z_sim(:, i) = cf.groupTransform(x_sim(:, i));
end

% --- Process Open Loop Data ---
N_sim_ol = length(t_sim_ol);
z_sim_ol = zeros(2*cf.n, N_sim_ol);
for i = 1:N_sim_ol
    z_sim_ol(:, i) = cf.groupTransform(x_sim_ol(:, i));
end

% --- Process Collocated Only Data ---
N_sim_col = length(t_sim_col);
z_sim_col = zeros(2*cf.n, N_sim_col);
for i = 1:N_sim_col
    z_sim_col(:, i) = cf.groupTransform(x_sim_col(:, i));
end

% ------------------------------
line_width = 3.0; 
line_style = "-";

% Show simulations
figure
% 1. Theta Actuated
subplot(4, 1, 1)
hold on
plot(t_sim_ol, z_sim_ol(1:cf.m, :), 'LineWidth', line_width, 'LineStyle', '-.', 'Color', red_ol)
plot(t_sim, z_sim(1:cf.m, :), 'LineWidth', line_width, 'LineStyle', line_style, 'Color', blue_sofft)
plot(t_sim_col, z_sim_col(1:cf.m, :), 'LineWidth', line_width, 'LineStyle', '--', 'Color', green_col)
hold off
grid on
xlabel("$t$ [s]", 'Interpreter', 'latex')
ylabel("$\theta_a$", 'Interpreter', 'latex')
set(gca, 'FontSize', font_size)
set(gca, 'GridLineWidth', grid_linewidth)
legend("Open-Loop", "Non-Collocated (Full)", "Collocated Only", 'Location', 'best')

% 2. Theta Unactuated
subplot(4, 1, 2)
hold on
plot(t_sim_ol, z_sim_ol(cf.m + 1:cf.n, :), 'LineWidth', line_width, 'LineStyle', '-.', 'Color', red_ol)
plot(t_sim, z_sim(cf.m + 1:cf.n, :), 'LineWidth', line_width, 'LineStyle', line_style, 'Color', blue_sofft)
plot(t_sim_col, z_sim_col(cf.m + 1:cf.n, :), 'LineWidth', line_width, 'LineStyle', '--', 'Color', green_col)
hold off
grid on
xlabel("$t$ [s]", 'Interpreter', 'latex')
ylabel("$\theta_u$", 'Interpreter', 'latex')
set(gca, 'FontSize', font_size)
set(gca, 'GridLineWidth', grid_linewidth)

% 3. Theta Dot Actuated
subplot(4, 1, 3)
hold on
plot(t_sim_ol, z_sim_ol(cf.n + 1: cf.n + cf.m, :), 'LineWidth', line_width, 'LineStyle', '-.', 'Color', red_ol)
plot(t_sim, z_sim(cf.n + 1: cf.n + cf.m, :), 'LineWidth', line_width, 'LineStyle', line_style, 'Color', blue_sofft)
plot(t_sim_col, z_sim_col(cf.n + 1: cf.n + cf.m, :), 'LineWidth', line_width, 'LineStyle', '--', 'Color', green_col)
hold off
grid on
xlabel("$t$ [s]", 'Interpreter', 'latex')
ylabel("$\dot{\theta}_a$", 'Interpreter', 'latex')
set(gca, 'FontSize', font_size)
set(gca, 'GridLineWidth', grid_linewidth)

% 4. Theta Dot Unactuated
subplot(4, 1, 4)
hold on
plot(t_sim_ol, z_sim_ol(cf.n + cf.m + 1:end, :), 'LineWidth', line_width, 'LineStyle', '-.', 'Color', red_ol)
plot(t_sim, z_sim(cf.n + cf.m + 1:end, :), 'LineWidth', line_width, 'LineStyle', line_style, 'Color', blue_sofft)
plot(t_sim_col, z_sim_col(cf.n + cf.m + 1:end, :), 'LineWidth', line_width, 'LineStyle', '--', 'Color', green_col)
hold off
grid on
xlabel("$t$ [s]", 'Interpreter', 'latex')
ylabel("$\dot{\theta}_u$", 'Interpreter', 'latex')
set(gca, 'FontSize', font_size)
set(gca, 'GridLineWidth', grid_linewidth)

%% Video
% T1.plotqt(t_sim, x_sim', "record", true, "video_name", "rsip_swingup_dist")

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
            % This function only uses Kpa and Kda inside
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
    % WrapToPi
    theta_tilde_a = wrapToPi(theta_des_a - theta_a);
    
    % Compute Control Law
    u = options.Kpa*theta_tilde_a - options.Kda*(theta_dot_a) + Ga + Ka;
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