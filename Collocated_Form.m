%%% Collocated Form for GVS Linkages %%%
classdef Collocated_Form < handle
    properties
        robot_linkage
    end

    methods
        function obj = Collocated_Form(robot_linkage)
            % Init
            obj.robot_linkage = robot_linkage;
        end

        % function theta = transform(q)
        %     % This function transform the coordinate q in the
        %     % collocated theta.
        % 
        %     % 
        % end
    end
end