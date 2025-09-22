%%% Collocated Form for GVS Linkages %%%
classdef Collocated_Form < handle
    properties
        robot_linkage
        n
        m
    end

    methods
        function obj = Collocated_Form(robot_linkage)
            % Init
            obj.robot_linkage = robot_linkage;
            obj.n = obj.robot_linkage.ndof;
            obj.m = obj.robot_linkage.nact;
        end

        %% Strain Field Function
        function xi = get_xi(obj, q, s)
            % Rescale
            L = obj.robot_linkage.VLinks(1).L;
            Phi_scale = diag([1/L 1/L 1/L 1 1 1]);
            
            % Compute strain
            xi = Phi_scale*obj.robot_linkage.CVRods{1}(2).Phi_h(s, obj.robot_linkage.CVRods{1}(2).Phi_dof, obj.robot_linkage.CVRods{1}(2).Phi_odr)*q + obj.robot_linkage.CVRods{1}(2).xi_starfn(s);
        end
        
        %% Actuation Matrix Function
        function [Btau, xi] = distributedActuationMatrix(obj, q, s)
            % Init
            Btau = zeros(6, obj.m);
        
            % Get xi
            xi = obj.get_xi(q, s);
            
            % Build Columns
            for i = 1:obj.m
                % Cable Routings
                d = obj.robot_linkage.CableFunction.dc_fn{i}(s);
                d_prime = obj.robot_linkage.CableFunction.dcp_fn{i}(s);
        
                % Tangent
                t = skew(xi(1:3))*d + xi(4:6) + d_prime;
                Btau(:, i) = [skew(d)*t; t];
            end
        end
        
        %% Actuators Length
        function La = actuatorLengths(obj, q)
            % Init
            La = zeros(obj.m, 1);
        
            % Gauss-Legendre Integration
            Xs = obj.robot_linkage.CVRods{1}(2).Xs;
            Ws = obj.robot_linkage.CVRods{1}(2).Ws;
        
            for i = 1:length(Xs)
                % Actuation Matrix
                [Btau, xi] = distributedActuationMatrix(obj, q, Xs(i));
                
                % Actuator Path
                for j = 1:obj.m
                    d_prime = obj.robot_linkage.CableFunction.dcp_fn{j}(Xs(i));
                    La(j) = La(j) + Ws(i)*(Btau(:, j)')*(xi + [zeros(3, 1); d_prime]);
                end
            end
        end

        %% Compute Transformation
        function theta = transform(obj, q)
            % This function transform the coordinate q in the
            % collocated theta.

            % For now only underactuated system
            assert(obj.m < obj.n, "Class supports only Underactuated system for now.");

            % Passive Output
            y = obj.actuatorLengths(q);
            
            % Change of Coordinates for Underactuated system
            theta = [y; zeros(obj.n - obj.m, 1)] + blkdiag(zeros(obj.m, obj.m), eye(obj.n - obj.m))*q;
        end

        %% Compute Transformation Dot
        function theta_dot = transformDot(obj, q, q_dot)
            theta_dot = (obj.robot_linkage.ActuationMatrix(q)')*q_dot;
        end
    end
end