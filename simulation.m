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
robot_name = "conical_hsupport";

% mat ext
file_name = "robot_linkage";
mat_ext = ".mat";

% Load Robot and Data
load(fullfile("robots", robot_name, "robot_linkage" + mat_ext));

%% Simulation Setup
fs = 1e+3;
tf = 5;
t0 = 0;
t = 0:(1/fs):tf;
N_sim = length(t);
n = T1.ndof;
q0 = randn(n, 1);
qdot0 = zeros(n, 1);
x0 = [q0; qdot0];

% Collocation Object
cf = Collocated_Form(T1);

% Desired Configuration
q_des = zeros(T1.ndof, 1);

% Store Result
x = zeros(2*n, N_sim);
u = zeros(T1.nact, N_sim);
x(:, 1) = x0;

% Visualization
x_sim = x0;
t_sim = t(1);

%% Main Cycle
for i = 1:(N_sim-1)
    % Compute Control Input
    u(:, i) = collocated_PD_FF(cf, q_des, x(1:n, i), x(n + 1:end, i));
    
    % Simulate
    ODEFUN = @(t, xk) dynamics(T1, t, xk, u(:, i));
    [TOUT, YOUT] = ode15s(ODEFUN, t(i:i+1), x(:, i));

    % Update State
    x(:, i + 1) = YOUT(end, :)';

    % Update Time and State | (Discard x(i) and t(i))
    t_sim = [t_sim, TOUT(2:end)'];
    x_sim = [x_sim, YOUT(2:end, :)'];

    % Iterations
    disp("Iteration " + num2str(i))
end

%% Plot Result
figure
subplot(2, 1, 1)
plot(t_sim, x_sim(1:n, :), 'LineWidth', 2.0)
hold on
yline(q_des, 'LineWidth', 1.5, 'LineStyle', '--')
hold off
grid on
xlabel("$t$ [s]", 'Interpreter', 'latex')
ylabel("$q$", 'Interpreter', 'latex')

subplot(2, 1, 2)
plot(t_sim, x_sim(n + 1:end, :), 'LineWidth', 2.0)
grid on
xlabel("$t$ [s]", 'Interpreter', 'latex')
ylabel("$\dot{q}$", 'Interpreter', 'latex')


function u = collocated_PD_FF(cf_obj, q_des, q, q_dot, Kpa, Kda)
    arguments
        % Collocation Form object: useful for conversion
        cf_obj

        % Desired state (in joint space, not in the collocated variables)
        q_des

        % Feedback terms
        q
        q_dot

        % Gains
        Kpa = ones(cf_obj.m);
        Kda = ones(cf_obj.m);
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
    u = Kpa*(theta_des_a - theta_a) - Kda*(theta_dot_a) + Ga + Ka;
end

function x_dot = dynamics(robot_linkage, t, x, u)
    % dynamicsSolver function
    [y,~,~,~] = robot_linkage.dynamicsSolver(t, x, u);
    
    % dxdt = [q_dot; q_2dot]
    x_dot = [x(robot_linkage.ndof + 1:end); y(1:robot_linkage.ndof)];
end