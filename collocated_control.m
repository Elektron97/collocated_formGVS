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

% Remove useless variables
% clear all

%% Load Data
robot_name = "conical_hsupport";

% mat ext
file_name = "robot_linkage";
mat_ext = ".mat";

% Load Robot and Data
load(fullfile("robots", robot_name, "robot_linkage" + mat_ext));

%% Collocated Form
cf = Collocated_Form(T1);

Btau = distributedActuationMatrix(T1, zeros(T1.ndof, 1), 0)

%% Actuation Matrix Function
function Btau = distributedActuationMatrix(robot_linkage, q, s)
    
    % Init
    m = robot_linkage.nact;
    Btau = zeros(6, m);
    L = robot_linkage.VLinks(1).L;
    Phi_scale = diag([1/L 1/L 1/L 1 1 1]);

    % Compute Strain
    xi = Phi_scale*robot_linkage.CVRods{1}(2).Phi_h(s, robot_linkage.CVRods{1}(2).Phi_dof, robot_linkage.CVRods{1}(2).Phi_odr)*q + robot_linkage.CVRods{1}(2).xi_starfn(s);
    
    % Build Columns
    for i = 1:m
        % Cable Routings
        d = robot_linkage.CableFunction.dc_fn{i}(s);
        d_prime = robot_linkage.CableFunction.dcp_fn{i}(s);

        % Tangent
        t = skew(xi(1:3))*d + xi(4:6) + d_prime;
        Btau(:, i) = [skew(d)*t; t];
    end
end
