%%% Collocated Form for GVS Linkages %%%
classdef Collocated_Form < handle
    properties
        robot_linkage
        n
        m
        L
        Phi_scale
    end

    methods
        function obj = Collocated_Form(robot_linkage)
            % Init
            obj.robot_linkage = robot_linkage;
            obj.n = obj.robot_linkage.ndof;
            obj.m = obj.robot_linkage.nact;

            % Useful Variables
            obj.L = obj.robot_linkage.VLinks(1).L;
            obj.Phi_scale =  diag([1/obj.L 1/obj.L 1/obj.L 1 1 1]);
        end

        %% Strain Field Function
        function xi = get_xi(obj, q, s)            
            % Compute strain
            xi = obj.Phi_scale*obj.robot_linkage.CVRods{1}(2).Phi_h(s, obj.robot_linkage.CVRods{1}(2).Phi_dof, obj.robot_linkage.CVRods{1}(2).Phi_odr)*q + obj.robot_linkage.CVRods{1}(2).xi_starfn(s);
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

        %% Actuation Matrix
        function A = actuationMatrix(obj, q)
            % Init
            A = zeros(obj.n, obj.m);
        
            % Gauss-Legendre Integration
            Xs = obj.robot_linkage.CVRods{1}(2).Xs;
            Ws = obj.robot_linkage.CVRods{1}(2).Ws;
        
            for i = 1:length(Xs)
                % Functional Basis
                Bq = obj.Phi_scale*obj.robot_linkage.CVRods{1}(2).Phi_h(Xs(i), obj.robot_linkage.CVRods{1}(2).Phi_dof, obj.robot_linkage.CVRods{1}(2).Phi_odr);

                % Actuation Matrix
                [Btau, ~] = distributedActuationMatrix(obj, q, Xs(i));
                
                % Actuator Path
                A = A + Ws(i)*(Bq')*Btau;
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
        function [theta, theta_dot] = transform(obj, q, qdot)
            % This function transform the coordinate q in the
            % collocated theta.

            % For now only underactuated system
            assert(obj.m <= obj.n, "Class supports only Underactuated system for now.");
            
            % Compute Actuation Matrix
            A = obj.actuationMatrix(q);
            assert(rank(A) == obj.m, "Actuation Matrix is not full rank.")

            % Passive Output
            y = obj.actuatorLengths(q);

            % To apply Pustina's Theorem, we need to reorder the rows of A
            % in such a way, the first m rows are full rank.
            % To perform this, we use QR algorithm.
            [~, ~, P] = qr(A');
            qp = (P')*q;
            qdotp = (P')*qdot;
            
            % Change of Coordinates for Underactuated system
            theta = [y; zeros(obj.n - obj.m, 1)] + blkdiag(zeros(obj.m, obj.m), eye(obj.n - obj.m))*qp;
            theta_dot = [(A')*qdot; [zeros(obj.n - obj.m, obj.m), eye(obj.n - obj.m)]*qdotp];
        end
    end
end