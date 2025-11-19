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

%% Gains
Kpa = 2.0e-2;
Kpu = 5.0e-2; % 5.0e-2

%% Plot Lyapunov
% Graphical Parameters
fontsize = 20;
grid_linewidth = 3;

% Colormap
load("SoFFTColormap.mat")

% Check if the robot has at least 2 DOFs for this plot
if n < 2
    fprintf('Robot has < 2 DOFs, skipping q1 vs q2 plot.\n');
else
    fprintf('Setting up grid for q1 vs q2 plot...\n');
    
    % --- MODIFIED SECTION START ---
    % Define the range strictly to +/- 1.5 around the equilibrium
    range_val = 2*pi;
    q1_range = linspace(q_des(1) - range_val, q_des(1) + range_val, 100); 
    q2_range = linspace(q_des(2) - range_val, q_des(2) + range_val, 100); 
    % --- MODIFIED SECTION END ---
    
    % 2. Create a 2D grid
    [Q1, Q2] = meshgrid(q1_range, q2_range);
    
    % 3. Initialize empty matrices
    V_grid = zeros(size(Q1));
    V_quad_grid = zeros(size(Q1)); 
    V_hessian_grid = zeros(size(Q1)); 
    
    % 4. Define the fixed values
    q_dot_fixed = q_dot_des; 
    
    % --- Computation ---
    [n_rows, n_cols] = size(Q1);
    
    fprintf('Calculating V values (this may take a moment)...\n');
    for i = 1:n_rows
        for j = 1:n_cols
            % Construct q vector
            q_val = q_des; 
            q_val(1) = Q1(i, j);
            q_val(2) = Q2(i, j);
            
            % Call function
            [V_grid(i, j), V_quad_grid(i, j), V_hessian_grid(i, j)] = ...
                lyapunov(cf, q_val, q_dot_fixed, q_des, q_dot_des, Kpa, Kpu);
        end
    end
    fprintf('Calculation complete.\n');
    
    % --- Visualization 1: Comparison (Full V vs Hessian) ---
    figure;
    hold on;
    
    % Plot 1: The Full Nonlinear Lyapunov Function (Surface)
    s1 = surf(Q1, Q2, V_grid);
    % s1.FaceAlpha = 0.8;        % Slightly transparent
    s1.EdgeColor = 'none';     % No lines for the smooth surface
    s1.DisplayName = 'Full Lyapunov V';
    shading interp;
    colormap(SoFFTColormap)
    colorbar;
    
    % % Plot 2: The Hessian Approximation (Mesh/Wireframe)
    % s2 = mesh(Q1, Q2, V_hessian_grid);
    % s2.FaceColor = 'none';     % See-through
    % % s2.EdgeColor = 'k';        % Black wireframe lines
    % s2.EdgeAlpha = 0.5;
    % s2.DisplayName = 'Hessian Approximation';

    % Mark the equilibrium
    [V_min, ~, ~] = lyapunov(cf, q_des, q_dot_des, q_des, q_dot_des, Kpa, Kpu);
    plot3(q_des(1), q_des(2), V_min, 'r*', 'MarkerSize', 15, 'LineWidth', 3, ...
        'DisplayName', 'Equilibrium');
    
    % Styling
    xlabel('$q_1$', 'Interpreter', 'latex');
    ylabel('$q_2$', 'Interpreter', 'latex');
    zlabel('Energy Value');
    title(['Comparison: Nonlinear V vs Hessian (Range \pm' num2str(range_val) ')']);
    legend('show', 'Location', 'best');
    view(3); 
    grid on;
    set(gca, 'FontSize', fontsize);
    set(gca, 'GridLineWidth', grid_linewidth);
    hold off;

    % --- Visualization 2: 2D Contour Plot (Total V) ---
    figure;
    contourf(Q1, Q2, V_grid, 20);
    xlabel('$q_1$', 'Interpreter', 'latex');
    ylabel('$q_2$', 'Interpreter', 'latex');
    title('Lyapunov Level Sets');
    colormap(SoFFTColormap)
    colorbar;
    axis equal;
    grid on;
    hold on;
    plot(q_des(1), q_des(2), 'r+', 'MarkerSize', 16, 'LineWidth', 3);
    legend('Level Sets', 'Equilibrium Point');
    hold off;
    set(gca, 'FontSize', fontsize)
    set(gca, 'GridLineWidth', grid_linewidth)
    
    % --- Visualization 3: 3D Surface (Quadratic Term ONLY) ---
    figure;
    surf(Q1, Q2, V_quad_grid);
    xlabel('$q_1$', 'Interpreter', 'latex');
    ylabel('$q_2$', 'Interpreter', 'latex');
    zlabel('Control Quadratic Term');
    colormap(SoFFTColormap)
    colorbar;
    shading interp;
    hold on;
    [~, V_quad_min, ~] = lyapunov(cf, q_des, q_dot_des, q_des, q_dot_des, Kpa, Kpu);
    plot3(q_des(1), q_des(2), V_quad_min, 'r*', 'MarkerSize', 15, 'LineWidth', 3);
    legend('Control Quadratic Term', 'Equilibrium Point');
    hold off;
    set(gca, 'FontSize', fontsize)
    set(gca, 'GridLineWidth', grid_linewidth)
end

%% Lyapunov Function Definition
function [V_total, V_quad, V_hessian] = lyapunov(cf_obj, q, q_dot, q_des, q_dot_des, Kpa, Kpu)
    % Potential Energy at the equilbrium
    [~, Uel_des, Ug_des] = cf_obj.mechanicalEnergy(q_des, q_dot_des);
    
    % Mechanical Energy at q, qdot
    [T, Uel, Ug] = cf_obj.mechanicalEnergy(q, q_dot);
    
    % Correction term
    [~, G_theta, K_theta, ~] = cf_obj.transformSystem(q_des, q_dot_des);
    [theta, ~] = cf_obj.transform(q, q_dot);
    
    % Desired Equilibrium in theta
    [theta_des, ~] = cf_obj.transform(q_des, q_dot_des);
    
    % Error in collocated variable
    theta_tilde = theta_des - theta;
    theta_tilde_a = theta_tilde(1:cf_obj.m);
    
    % Gain
    Kp = [Kpa, Kpu; zeros(cf_obj.p, cf_obj.n)];
    
    % --- Break down the V terms ---
    % 1. Mechanical Energy part
    V_energy = T + Uel + Ug - Uel_des - Ug_des;
    % 2. Correction part
    V_correction = (theta_tilde_a')*(G_theta(1:cf_obj.m) + K_theta(1:cf_obj.m));
    % 3. Quadratic part
    V_quad = 0.5*(theta_tilde')*Kp*theta_tilde;
    
    % Compute total Lyapunov Function
    V_total = V_energy + V_correction + V_quad;
    
    % --- Potential Hessian --- %
    % Note: Ensure potential_hessian is available in your path or defined
    Gamma = potential_hessian(cf_obj.robot_linkage, q_des, zeros(cf_obj.n, 1), zeros(cf_obj.m, 1));
    Jh = cf_obj.jacobian(q);
    
    % Check for singularity or inversion issues if needed
    % Gamma_theta = (inv(Jh)')*Gamma*inv(Jh); 
    % Better numerical stability than inv(Jh):
    Gamma_theta = (Jh') \ (Gamma / Jh); 
    
    % Quadratic Form for Hessian approximation
    V_hessian = 0.5*(theta_tilde')*(Gamma_theta + Kp)*theta_tilde;
end