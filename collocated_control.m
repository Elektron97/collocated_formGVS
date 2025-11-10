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

% mat ext
file_name = "robot_linkage";
mat_ext = ".mat";

% Load Robot and Data
load(fullfile("robots", robot_name, "robot_linkage" + mat_ext));

%% Collocated Form
cf = Collocated_Form(T1);

%% Equilibria of the System
% Load EquilibriaGVS repo
addpath(fullfile("..", "GVS-OptimalControl", "EquilibriaGVS"))
addpath(fullfile("..", "GVS-OptimalControl", "EquilibriaGVS", "functions"))

% Call equilibriaGVS function
u = zeros(T1.nact, 1);
q = equilibriaGVS(T1, "input", u);

% WrapToPi correction for the revolut joint
if T1.CVRods{1}(1).dof >= 1
    q(1, :) = wrapToPi(q(1, :));
end

% Filter Equilibria
q = filterEquilibria(q);

%% Lyapunov Function
% function Vc = collocated_Lyapunov(theta_a, theta_ad, Ga_eq, Ka_eq, K_pa)
%     arguments
%         theta_a
%         theta_ad
%         Ga_eq
%         Ka_eq
%         K_pa = eye(length(theta_a));
%     end
% 
%     % Compute Collocated Lyapunov Term
%     theta_tilde_a = theta_ad - theta_a;
% 
%     % Minimum Eigenvalue from the matrix
%     [lambda, ~] = eig(K_pa);
%     lambda_min = min(lambda);
% 
%     Vc = 0.5*(theta_tilde_a')*K_pa*theta_tilde_a + (theta_tilde_a')*(Ga_eq + Ka_eq) + 0.5*(norm(Ga_eq + Ka_eq)^2)/lambda_min;
% end

% function H = energy(cf_obj, theta, theta_dot)
%     arguments
%         cf_obj,
%         theta,
%         theta_ad
%         Ga_eq
%         Ka_eq
%         K_pa = eye(length(theta_a));
%     end
% 
%     % Compute Collocated Lyapunov Term
%     theta_tilde_a = theta_ad - theta_a;
% 
%     % Minimum Eigenvalue from the matrix
%     [lambda, ~] = eig(K_pa);
%     lambda_min = min(lambda);
% 
%     Vc = 0.5*(theta_tilde_a')*K_pa*theta_tilde_a + (theta_tilde_a')*(Ga_eq + Ka_eq) + 0.5*(norm(Ga_eq + Ka_eq)^2)/lambda_min;
% end