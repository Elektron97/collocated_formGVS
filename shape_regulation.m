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
% 
% T1 = T1.Update();
% 
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

if is_stable
    set(gca, 'XScale', 'log')
end

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

%% Apply LQR
% qa = 1e+3;
% qu =  1e+0;
% qad = 1e+3;
% qud = 1e+0;
% r = 1e+0; % 1e-3

% Brayson Rule
qa = (1/2)^2;
qu = (1/3)^2;
qad = (1/2)^2;
qud = (1/10)^2;
r = 1e-3; % 1e-3

% Build Weight Matrices
Q = blkdiag(qa*eye(cf.m), qu*eye(cf.p), qad*eye(cf.m), qud*eye(cf.p));
R = r*eye(cf.m);

% Solve Riccati
[K_lqr, P, lambda_cl] = lqr(A_theta, B_theta, Q, R);

% Extract Gains
Kpa = K_lqr(1:cf.m);
Kpu = K_lqr(cf.m + 1:cf.n);
Kda = K_lqr(cf.n + 1: cf.n + cf.m);
Kdu = K_lqr(cf.n + cf.m + 1:end);

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
% ---------------------------------------------------------
% 1. Simulate Non-Collocated (Full Gains: Kpa, Kda, Kpu, Kdu)
control_law = "noncollocated_PD_FF";
ODEFUN_FULL = @(t, xk) closed_loop(cf, t, xk, "q_des", q_des, ...
                                "control_law", control_law, ...
                                "Kpa", Kpa, "Kda", Kda, ...
                                "Kpu", Kpu, "Kdu", Kdu);
[t_sim, x_sim] = ode15s(ODEFUN_FULL, t, x0);
x_sim = x_sim';

% ---------------------------------------------------------
% 2. Simulate Collocated Only (Partial Gains: Kpa, Kda only)
control_law_col = "collocated_PD_FF";
ODEFUN_COL = @(t, xk) closed_loop(cf, t, xk, "q_des", q_des, ...
                                "control_law", control_law_col, ...
                                "Kpa", Kpa, "Kda", Kda);
[t_sim_col, x_sim_col] = ode15s(ODEFUN_COL, t, x0);
x_sim_col = x_sim_col';

% ---------------------------------------------------------
% 3. Simulate Non-Collocated with Kpu = 0
control_law_nokpu = "damping_injection";
ODEFUN_NO_KPU = @(t, xk) closed_loop(cf, t, xk, "q_des", q_des, ...
                                "control_law", control_law_nokpu, ...
                                "Kpa", Kpa, "Kda", Kda); 
[t_sim_nokpu, x_sim_nokpu] = ode15s(ODEFUN_NO_KPU, t, x0);
x_sim_nokpu = x_sim_nokpu';

%% Visualization
% Collocated Desired variables
z_des = cf.groupTransform([q_des; zeros(cf.n, 1)]);

% --- Process Non-Collocated Data (Full) ---
N_sim = length(t_sim);
z_sim = zeros(2*cf.n, N_sim);
for i = 1:N_sim
    z_sim(:, i) = cf.groupTransform(x_sim(:, i));
end
% --- Process Collocated Only Data ---
N_sim_col = length(t_sim_col);
z_sim_col = zeros(2*cf.n, N_sim_col);
for i = 1:N_sim_col
    z_sim_col(:, i) = cf.groupTransform(x_sim_col(:, i));
end
% --- Process Non-Collocated (No Kpu) Data ---
N_sim_nokpu = length(t_sim_nokpu);
z_sim_nokpu = zeros(2*cf.n, N_sim_nokpu);
for i = 1:N_sim_nokpu
    z_sim_nokpu(:, i) = cf.groupTransform(x_sim_nokpu(:, i));
end
% ------------------------------
line_width = 4.0; 
line_style = "-";
% Show simulations
figure
tiledlayout(2, 2, 'TileSpacing','Compact', 'Padding', 'Compact')

% 1. Theta Actuated
% subplot(2, 2, 1)
nexttile
hold on
plot(t_sim, z_sim(1:cf.m, :), 'LineWidth', line_width, 'LineStyle', line_style, 'Color', blue_sofft)
if ~T1.CEF
    plot(t_sim_col, z_sim_col(1:cf.m, :), 'LineWidth', line_width, 'LineStyle', '-.', 'Color', green_col)
    % plot(t_sim_nokpu, z_sim_nokpu(1:cf.m, :), 'LineWidth', line_width, 'LineStyle', '-.', 'Color', yellow_nokpu)
end
yline(z_des(1:cf.m), 'LineWidth', 3.0, 'LineStyle', '--', 'Color', grey_mid)
hold off
grid on
xlim([t_sim(1), t_sim(end)])
xlabel("$t$ [s]", 'Interpreter', 'latex')
ylabel("$\theta_a$", 'Interpreter', 'latex')
set(gca, 'FontSize', font_size)
set(gca, 'GridLineWidth', grid_linewidth)
% legend("Non-Collocated", "Collocated Only", "Damping Injection", 'Interpreter', 'latex', 'Location', 'northeast')
% legend("Non-Collocated PD+FF", "Collocated PD+FF", "Target", 'Interpreter', 'latex', 'Location', 'best')

% 2. Theta Unactuated
% subplot(2, 2, 2)
nexttile
hold on
plot(t_sim, z_sim(cf.m + 1:cf.n, :), 'LineWidth', line_width, 'LineStyle', line_style, 'Color', blue_sofft)
if ~T1.CEF
    plot(t_sim_col, z_sim_col(cf.m + 1:cf.n, :), 'LineWidth', line_width, 'LineStyle', '-.', 'Color', green_col)
    % plot(t_sim_nokpu, z_sim_nokpu(cf.m + 1:cf.n, :), 'LineWidth', line_width, 'LineStyle', '-.', 'Color', yellow_nokpu)
end
yline(z_des(cf.m + 1:cf.n), 'LineWidth', 3.0, 'LineStyle', '--', 'Color', grey_mid)
hold off
grid on
xlim([t_sim(1), t_sim(end)])
xlabel("$t$ [s]", 'Interpreter', 'latex')
ylabel("$\theta_u$", 'Interpreter', 'latex')
set(gca, 'FontSize', font_size)
set(gca, 'GridLineWidth', grid_linewidth)
% legend("Non-Collocated PD+FF", "Collocated PD+FF", 'Interpreter', 'latex', 'Location', 'northeast')

% 3. Theta Dot Actuated
% subplot(2, 2, 3)
nexttile
hold on
plot(t_sim, z_sim(cf.n + 1: cf.n + cf.m, :), 'LineWidth', line_width, 'LineStyle', line_style, 'Color', blue_sofft)
if ~T1.CEF
    plot(t_sim_col, z_sim_col(cf.n + 1: cf.n + cf.m, :), 'LineWidth', line_width, 'LineStyle', '-.', 'Color', green_col)
    % plot(t_sim_nokpu, z_sim_nokpu(cf.n + 1: cf.n + cf.m, :), 'LineWidth', line_width, 'LineStyle', '-.', 'Color', yellow_nokpu)
end
yline(z_des(cf.n + 1: cf.n + cf.m), 'LineWidth', 3.0, 'LineStyle', '--', 'Color', grey_mid)
hold off
grid on
xlim([t_sim(1), t_sim(end)])
xlabel("$t$ [s]", 'Interpreter', 'latex')
ylabel("$\dot{\theta}_a$", 'Interpreter', 'latex')
set(gca, 'FontSize', font_size)
set(gca, 'GridLineWidth', grid_linewidth)
% legend("Non-Collocated PD+FF", "Collocated PD+FF", 'Interpreter', 'latex', 'Location', 'northeast')

% 4. Theta Dot Unactuated
% subplot(2, 2, 4)
nexttile
hold on
plot(t_sim, z_sim(cf.n + cf.m + 1:end, :), 'LineWidth', line_width, 'LineStyle', line_style, 'Color', blue_sofft)
if ~T1.CEF
    plot(t_sim_col, z_sim_col(cf.n + cf.m + 1:end, :), 'LineWidth', line_width, 'LineStyle', '-.', 'Color', green_col)
    % plot(t_sim_nokpu, z_sim_nokpu(cf.n + cf.m + 1:end, :), 'LineWidth', line_width, 'LineStyle', '-.', 'Color', yellow_nokpu)
end
yline(z_des(cf.n + cf.m + 1), 'LineWidth', 3.0, 'LineStyle', '--', 'Color', grey_mid)
hold off
grid on
xlim([t_sim(1), t_sim(end)])
xlabel("$t$ [s]", 'Interpreter', 'latex')
ylabel("$\dot{\theta}_u$", 'Interpreter', 'latex')
set(gca, 'FontSize', font_size)
set(gca, 'GridLineWidth', grid_linewidth)
% legend("Non-Collocated PD+FF", "Collocated PD+FF", 'Interpreter', 'latex', 'Location', 'northeast')

%% Reconstruct the Input (Debug)
u_sim = zeros(cf.m, length(t_sim));

for i = 1:length(t_sim)
    u_sim(:, i) = noncollocated_PD_FF(cf, q_des, x_sim(1:cf.n, i), x_sim(cf.n + 1:end, i), ...
                            "Kpa", Kpa, "Kda", Kda, ...
                            "Kpu", Kpu, "Kdu", Kdu);
end

figure
plot(t_sim, u_sim, 'LineWidth', line_width, 'Color', blue_sofft)
grid on
xlabel("$t$ [s]")
ylabel("$u$ [N $\cdot$ m]")

%% Save Figure
figure_path = fullfile("figures", robot_name);
save_figure = false;

if save_figure
    file_name = "swingup_simulation_dist";
    % savefig(fullfile(figure_path, file_name + ".fig"));

    % Rendering SVG
    set(gcf,'renderer','painters');
    saveas(gcf, fullfile(figure_path, file_name + ".svg"));
end

%% Snapshots
% % Mid Points
% T1.PlotParameters.ClosePrevious = true;
% T1.VLinks.alpha = 0.2;
% % T1.VLinks.color = hex2rgb(grey_mid);
% 
% if ~T1.CEF
%     frames = [1, 50, 100:280:1000, length(t_sim)];
% else
%     % frames = [4800:100:5000, 5000:100:5500, 6000];
%     frames = 4800:250:6200;
% end
% % Colors
% colors = SoFFTColormap(length(frames));
% 
% figure
% plot(nan, nan)
% hold on
% for i = 2:length(frames)-1
%     T1.VLinks.color = colors(i,:);
%     T1.plotq(x_sim(1:T1.ndof, frames(i))');
%     T1.PlotParameters.ClosePrevious = false;
% end
% hold off
% set(gca, 'FontSize', 25)
% 
% if save_figure
%     file_name = "swingup_frames_dist";
%     % savefig(fullfile(figure_path, file_name + ".fig"));
% 
%     % Rendering SVG
%     set(gcf,'renderer','painters');
%     saveas(gcf, fullfile(figure_path, file_name + ".svg"));
% end
% 
% figure
% plot(nan)
% hold on
% T1.PlotParameters.ClosePrevious = true;
% T1.VLinks.alpha = 1.0;
% T1.VLinks.color = colors(2, :);
% T1.plotq(x_sim(1:T1.ndof, frames(3))');
% 
% T1.PlotParameters.ClosePrevious = false;
% % 
% % T1.VLinks.color = colors(end, :);
% % T1.plotq(x_sim(1:T1.ndof, frames(end))');
% hold off
% set(gca, 'FontSize', 25)
% 
% if save_figure
%     file_name = "swingup_target_dist";
%     % savefig(fullfile(figure_path, file_name + ".fig"));
% 
%     % Rendering SVG
%     set(gcf,'renderer','painters');
%     saveas(gcf, fullfile(figure_path, file_name + ".svg"));
% end
% 
% if ~T1.CEF
%     % Mid Points
%     T1.PlotParameters.ClosePrevious = true;
%     T1.VLinks.alpha = 0.2;
%     % T1.VLinks.color = hex2rgb(grey_mid);
% 
%     frames = [1, 50, 100:280:1000, length(t_sim_col)];
%     % Colors
%     colors = SoFFTColormap(length(frames));
% 
%     figure
%     plot(nan, nan)
%     hold on
%     for i = 2:length(frames)-1
%         T1.VLinks.color = colors(i,:);
%         T1.plotq(x_sim_col(1:T1.ndof, frames(i))');
%         T1.PlotParameters.ClosePrevious = false;
%     end
%     hold off
%     set(gca, 'FontSize', 25)
% 
%     if save_figure
%         file_name = "swingup_frames_col";
%         % savefig(fullfile(figure_path, file_name + ".fig"));
% 
%         % Rendering SVG
%         set(gcf,'renderer','painters');
%         saveas(gcf, fullfile(figure_path, file_name + ".svg"));
%     end
% 
%     figure
%     plot(nan)
%     hold on
%     T1.PlotParameters.ClosePrevious = true;
%     T1.VLinks.alpha = 1.0;
%     T1.VLinks.color = hex2rgb(blue_sofft);
%     T1.plotq(x_sim_col(1:T1.ndof, 1)');
% 
%     T1.PlotParameters.ClosePrevious = false;
% 
%     T1.VLinks.color = hex2rgb(yellow_nokpu);
%     T1.plotq(x_sim_col(1:T1.ndof, end)');
%     hold off
%     set(gca, 'FontSize', 25)
% 
%     if save_figure
%         file_name = "swingup_target_col";
%         % savefig(fullfile(figure_path, file_name + ".fig"));
% 
%         % Rendering SVG
%         set(gcf,'renderer','painters');
%         saveas(gcf, fullfile(figure_path, file_name + ".svg"));
%     end
% end

%% Video
% T1.PlotParameters.ClosePrevious = true;
% T1.VLinks.alpha = 1.0;
% T1.VLinks.color = hex2rgb(blue_sofft);
% T1.plotqt(t_sim, x_sim', "record", true, "video_name", robot_name)

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

        case "damping_injection"
            u = damping_injection(cf, options.q_des, x(1:cf.n), x(cf.n + 1:end), "Kpa", options.Kpa, "Kda", options.Kda);

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
    % theta_tilde_a = wrapToPi(theta_des_a - theta_a);
    theta_tilde_a = theta_des_a - theta_a;
    
    % Compute Control Law
    u = options.Kpa*theta_tilde_a - options.Kda*(theta_dot_a) + Ga + Ka;
end

%% Non-collocated Law
function u = damping_injection(cf_obj, q_des, q, q_dot, options)
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
    [u_c, ~, ~, theta_dot] = collocated_PD_FF(cf_obj, q_des, q, q_dot, "Kpa", options.Kpa, "Kda", options.Kda);
    
    % Compute (optimal) Damping Injection
    [~, ~, ~, D_theta] = cf_obj.transformSystem(q_des, zeros(cf_obj.n, 1));
    Kdu = -2.0*D_theta(1:cf_obj.m, (cf_obj.m+1):end);


    % Kdu = -2 Dau
    u_nc = - Kdu*theta_dot(cf_obj.m + 1:end);
    
    % Compose the Actions
    u = u_c + u_nc;
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
        options.Kpu = zeros(cf_obj.m, cf_obj.p);
        options.Kdu = zeros(cf_obj.m, cf_obj.p);
    end
    
    % Collocated Term
    [u_c, theta_des, theta, theta_dot] = collocated_PD_FF(cf_obj, q_des, q, q_dot, "Kpa", options.Kpa, "Kda", options.Kda);
    
    % Add noncollocated term
    theta_tilde_u = theta_des(cf_obj.m + 1:end) - theta(cf_obj.m + 1:end);
    u_nc = options.Kpu*theta_tilde_u - options.Kdu*theta_dot(cf_obj.m + 1:end);
    
    % Compose the Actions
    u = u_c + u_nc;
end

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