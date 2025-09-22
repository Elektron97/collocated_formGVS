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

% Test change of coordinates
q0 = zeros(T1.ndof, 1);
theta = cf.transform(q0);