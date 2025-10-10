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
% robot_name = "rsip";
robot_name = "conical_hsupport";

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
T1.PlotParameters.XLim = [-0.6, 0.6];
T1.PlotParameters.YLim = [-0.6, 0.6];

% Colors
blue_sofft = "#086788";
red_target = "#f06543";
grey_mid = "#858583";
T1.VLinks.color = hex2rgb(blue_sofft);

% Camera Position
T1.PlotParameters.CameraPosition = [-0.0316   -0.0001   -6.9004];

% Damping Joint
T1.D(1, 1) = 1e-2;

% Update Linkage
T1 = T1.Update();

%% Simulation Setup
fs = 1e+3;
tf = 50;
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

% Method
control_law = "noncollocated_PD_FF";
% control_law = "collocated_PD_FF";
% control_law = "classic_PD_FF";

%% Simulate
ODEFUN = @(t, xk) closed_loop(cf, t, xk, "q_des", q_des, "control_law", control_law, ...
                                    "Kpa", 20*eye(cf.m),         "Kda", 10*eye(cf.m), ...
                                    "Kpu", 20*eye(cf.n - cf.m),  "Kdu", 10*eye(cf.n - cf.m));
[t_sim, x_sim] = ode15s(ODEFUN, t, x0);

% Column Notation
x_sim = x_sim';

%% Video
T1.plotqt(t_sim, x_sim', "record", false);

%% Plot Result (Configuration Space)
figure
subplot(4, 1, 1)
plot(t_sim, x_sim(1, :), 'LineWidth', 2.0)
hold on
yline(q_des(1), 'LineWidth', 1.5, 'LineStyle', '--')
hold off
grid on
xlabel("$t$ [s]", 'Interpreter', 'latex')
ylabel("$\theta$", 'Interpreter', 'latex')

subplot(4, 1, 2)
plot(t_sim, x_sim(2:end, :), 'LineWidth', 2.0)
hold on
yline(q_des(2:end), 'LineWidth', 1.5, 'LineStyle', '--')
hold off
grid on
xlabel("$t$ [s]", 'Interpreter', 'latex')
ylabel("$\kappa$", 'Interpreter', 'latex')

subplot(4, 1, 3)
plot(t_sim, x_sim(n + 1, :), 'LineWidth', 2.0)
grid on
xlabel("$t$ [s]", 'Interpreter', 'latex')
ylabel("$\dot{\theta}$", 'Interpreter', 'latex')

subplot(4, 1, 4)
plot(t_sim, x_sim(n + 2:end, :), 'LineWidth', 2.0)
grid on
xlabel("$t$ [s]", 'Interpreter', 'latex')
ylabel("$\dot{\kappa}$", 'Interpreter', 'latex')

function u = classic_PD_FF(cf_obj, q_des, q, q_dot, options)
    arguments
        % Collocation Form object: useful for conversion
        cf_obj

        % Desired state (in joint space, not in the collocated variables)
        q_des

        % Feedback terms
        q
        q_dot

        % Gains
        options.alpha = eye(cf_obj.m);
        options.beta = eye(cf_obj.m);
    end

    % Gravity and Stiffness Components
    [~, Geq, Keq, ~] = cf_obj.dynamicMarices(q_des, zeros(cf_obj.n, 1));

    % Compute Control Law
    A = cf_obj.actuationMatrix(q);
    u = pinv(A)*(Geq + Keq) + (options.alpha*A')*(q_des - q) - (options.beta*A')*q_dot;
end

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

    %% Additional Informations
    if nargout > 1
        
    end
end

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
        options.Kpu = eye(cf_obj.n - cf_obj.m);
        options.Kdu = eye(cf_obj.n - cf_obj.m);
    end

    % Underactuated Num. of variables
    p = cf_obj.n - cf_obj.m;

    % Collocated Term
    [u_c, theta_des, theta, theta_dot] = collocated_PD_FF(cf_obj, q_des, q, q_dot, "Kpa", options.Kpa, "Kda", options.Kda);

    % Add noncollocated term
    S = [zeros(cf_obj.m, p); eye(p)];
    Apinv = pinv(cf_obj.actuationMatrix(q));
    theta_tilde_u = theta_des(cf_obj.m + 1:end) - theta(cf_obj.m + 1:end);
    u_nc = Apinv*S*( options.Kpu*theta_tilde_u - options.Kdu*theta_dot(cf_obj.m + 1:end));

    % Compose the Actions
    u = u_c + u_nc;
end

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
        options.control_law = "classic_PD_FF";
        options.Kpa = eye(cf.robot_linkage.nact);
        options.Kda = eye(cf.robot_linkage.nact);
        options.Kpu = eye(cf.robot_linkage.ndof - cf.robot_linkage.nact);
        options.Kdu = eye(cf.robot_linkage.ndof - cf.robot_linkage.nact);
        options.q_des = zeros(2*cf.robot_linkage.ndof, 1);
    end

    % Compute Control Action
    u = zeros(cf.robot_linkage.nact, 1);

    % Select control law
    switch options.control_law
        case "classic_PD_FF"
            u = classic_PD_FF(cf, options.q_des, x(1:cf.n), x(cf.n + 1:end), "alpha", options.Kpa, "beta", options.Kda);
        case "collocated_PD_FF"
            u = collocated_PD_FF(cf, options.q_des, x(1:cf.n), x(cf.n + 1:end), "Kpa", options.Kpa, "Kda", options.Kda);
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