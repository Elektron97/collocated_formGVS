%%% Build Robot %%%
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

%% Useful Paths
% robot_name =  "rsip";
robot_name =  "rsip_extreme";
% robot_name =  "free_rod";
robot_dir = fullfile("robots", robot_name);

% Create Directory
if ~exist(robot_dir, 'dir')
   mkdir(robot_dir)
end

%% Flags
rebuild_link = true;
rebuild_linkage = true;

% File Names
mat_ext = ".mat";
link_file = "robot_link";
linkage_file = "robot_linkage";

%% Create or Load SoRoSimLink
if rebuild_link
    L1 = SorosimLink();
    save(fullfile(robot_dir, link_file + mat_ext), "L1");
else
    if exist(fullfile(robot_dir, link_file + mat_ext), 'file')
        load(fullfile(robot_dir, link_file + mat_ext));
    else
        disp("Error: Linkage not created.");
    end
end

%% Create or Load SoRoSimLinkage
if rebuild_linkage
    T1 = SorosimLinkage(L1);
    save(fullfile(robot_dir, linkage_file + mat_ext), "T1");
else
    if exist(fullfile(robot_dir, linkage_file + mat_ext), 'file')
        load(fullfile(robot_dir, linkage_file + mat_ext));
    else
        disp("Error: Linkage not created.");
    end
end

%% Improve Visualization
% Camera Position
T1.PlotParameters.Light = false;
T1.PlotParameters.ClosePrevious = false;

% Axes Limits
T1.PlotParameters.XLim = [-T1.VLinks.L, T1.VLinks.L];
T1.PlotParameters.YLim = [-T1.VLinks.L, T1.VLinks.L];

% Colors
blue_sofft = "#086788";
T1.VLinks.color = hex2rgb(blue_sofft);

% Camera Position
T1.PlotParameters.CameraPosition = [-0.0316   -0.0001   -6.9004];

% Update Linkage
T1 = T1.Update();

% Damping JointT1.D
T1.D(1, 1) = 0.5;

% Save Robot after modification
 save(fullfile(robot_dir, linkage_file + mat_ext), "T1");

%% Show Robot
% T1.plotq(0.2*randn(T1.ndof, 1))