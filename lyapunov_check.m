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
Kpa = 1;
Kpu = 2;

%% Plot Lyapunov
% We will visualize V by sweeping q(1) and q(2), while holding all
% other q_i and ALL q_dot at their desired equilibrium values.
% Check if the robot has at least 2 DOFs for this plot
if n < 2
    fprintf('Robot has < 2 DOFs, skipping q1 vs q2 plot.\n');
    % You might want to 'return' or just skip this section
else
    fprintf('Setting up grid for q1 vs q2 plot...\n');
    
    % 1. Define a reasonable range for q1 and q2
    % (Centered around the equilibrium q_des)
    % You may need to adjust the '0.5' to a range that suits your robot
    q1_range = linspace(q_des(1) -2*pi, q_des(1) + 2*pi, 100); % 100 points
    q2_range = linspace(q_des(2) - 10, q_des(2) + 10, 100); % 100 points
    
    % 2. Create a 2D grid
    [Q1, Q2] = meshgrid(q1_range, q2_range);
    
    % 3. Initialize empty matrices to store V and its quadratic part
    V_grid = zeros(size(Q1));
    V_quad_grid = zeros(size(Q1)); % <-- ADDED
    
    % 4. Define the fixed values (all velocities are zero)
    q_dot_fixed = q_dot_des; 
    
    % --- Computation ---
    [n_rows, n_cols] = size(Q1);
    
    fprintf('Calculating V values (this may take a moment)...\n');
    for i = 1:n_rows
        for j = 1:n_cols
            % Start with the equilibrium q vector
            q_val = q_des; 
            
            % Overwrite the first two elements with grid values
            q_val(1) = Q1(i, j);
            q_val(2) = Q2(i, j);
            
            % q_dot_val is held constant at q_dot_des (zeros)
            
            % Call your non-vectorizable function
            % <-- MODIFIED to get two outputs
            [V_grid(i, j), V_quad_grid(i, j)] = lyapunov(cf, q_val, q_dot_fixed, q_des, q_dot_des, Kpa, Kpu);
        end
    end
    fprintf('Calculation complete.\n');
    
    % --- Visualization 1: 3D Surface (Total V) ---
    figure;
    surf(Q1, Q2, V_grid);
    xlabel('Position (q_1)');
    ylabel('Position (q_2)');
    zlabel('Lyapunov Function Value (V)');
    title('Lyapunov Function Slice (VelocITIES = 0)');
    colorbar;
    shading interp;
    hold on;
    % Mark the equilibrium
    [V_min, ~] = lyapunov(cf, q_des, q_dot_des, q_des, q_dot_des, Kpa, Kpu);
    plot3(q_des(1), q_des(2), V_min, 'r*', 'MarkerSize', 10, 'LineWidth', 2);
    legend('Full Lyapunov V', 'Equilibrium Point');
    hold off;
    
    % --- Visualization 2: 2D Contour Plot (Total V) ---
    figure;
    contourf(Q1, Q2, V_grid, 20);
    xlabel('Position (q_1)');
    ylabel('Position (q_2)');
    title('Lyapunov Level Sets (Velocities = 0)');
    colorbar;
    axis equal;
    grid on;
    hold on;
    % Mark the equilibrium
    plot(q_des(1), q_des(2), 'r+', 'MarkerSize', 12, 'LineWidth', 2);
    legend('Level Sets', 'Equilibrium Point');
    hold off;
    
    % --- Visualization 3: 3D Surface (Quadratic Term ONLY) ---
    % <-- NEW PLOT ADDED
    figure;
    surf(Q1, Q2, V_quad_grid);
    xlabel('Position (q_1)');
    ylabel('Position (q_2)');
    zlabel('Lyapunov Value (V_{quad})');
    title('Quadratic Term (0.5*\theta^T*Kp*\theta) Only');
    colorbar;
    shading interp;
    hold on;
    % Mark the equilibrium
    [~, V_quad_min] = lyapunov(cf, q_des, q_dot_des, q_des, q_dot_des, Kpa, Kpu);
    plot3(q_des(1), q_des(2), V_quad_min, 'r*', 'MarkerSize', 10, 'LineWidth', 2);
    legend('Quadratic Term', 'Equilibrium Point');
    hold off;

end

%% Lyapunov Function Definition
% <-- MODIFIED to return two values
function [V_total, V_quad] = lyapunov(cf_obj, q, q_dot, q_des, q_dot_des, Kpa, Kpu)
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
    % 3. Quadratic part (the one you want to see separately)
    V_quad = 0.5*(theta_tilde')*Kp*theta_tilde;
    
    % Compute total Lyapunov Function
    V_total = V_energy + V_correction + V_quad;
end