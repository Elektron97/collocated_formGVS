%%% Collocated Form for GVS Linkages %%%
classdef Collocated_Form < handle
    properties
        robot_linkage
        n
        m
        joint_dof
        strain_dof
        n_sact
        L
        Phi_scale
        q_sing
    end

    methods
        function obj = Collocated_Form(robot_linkage)
            % Init
            obj.robot_linkage = robot_linkage;
            obj.n = obj.robot_linkage.ndof;
            obj.m = obj.robot_linkage.nact;

            % Joint and Soft DoFs
            obj.joint_dof = obj.robot_linkage.CVRods{1}(1).dof;
            obj.strain_dof = obj.robot_linkage.CVRods{1}(2).dof;

            % Soft Actuators
            obj.n_sact = obj.robot_linkage.n_sact;

            % Useful Variables
            obj.L = obj.robot_linkage.VLinks(1).L;
            obj.Phi_scale =  diag([1/obj.L 1/obj.L 1/obj.L 1 1 1]);

            % Test on Actuation Matrix
            obj.testActuation();
        end

        %% Strain Field Function
        function xi = get_xi(obj, q, s)
            % Assert for current implementation
            assert(isscalar(obj.robot_linkage.CVRods), "For now only single link linkage is supproted. Sorry, I'm lazy.")

            % Extract Joint and Strain q
            rod = obj.robot_linkage.CVRods{1};
            obj.joint_dof = rod(1).dof;
            obj.strain_dof = rod(2).dof;

            % Compute strain
            xi = obj.Phi_scale*rod(2).Phi_h(s, rod(2).Phi_dof, rod(2).Phi_odr)*q(obj.joint_dof + 1: obj.joint_dof + obj.strain_dof) + obj.robot_linkage.CVRods{1}(2).xi_starfn(s);
        end
        
        %% Distributed Actuation Matrix Function
        function [Btau, xi] = distributedActuationMatrix(obj, q, s)
            % Init
            Btau = zeros(6, obj.robot_linkage.n_sact);
        
            % Get xi
            xi = obj.get_xi(q, s);
            
            % Build Columns
            for i = 1:obj.robot_linkage.n_sact
                % Cable Routings
                d = obj.robot_linkage.CableFunction.dc_fn{i}(s);
                d_prime = obj.robot_linkage.CableFunction.dcp_fn{i}(s);
        
                % Tangent
                t = skew(xi(1:3))*d + xi(4:6) + d_prime;
                t = t / norm(t);
                Btau(:, i) = [skew(d)*t; t];
            end
        end

        %% Actuation Matrix
        function [A, P, Aa, Au] = actuationMatrix(obj, q)
            % Init
            A = zeros(obj.n, obj.m);

            % Joint Actuation Matrix
            if obj.joint_dof > 0
                A(:, 1) = obj.robot_linkage.Bj1;
            end
        
            if obj.n_sact > 0
                % Gauss-Legendre Integration
                Xs = obj.robot_linkage.CVRods{1}(2).Xs;
                Ws = obj.robot_linkage.CVRods{1}(2).Ws;
            
                for i = 1:length(Xs)
                    % Functional Basis
                    Bq = obj.Phi_scale*obj.robot_linkage.CVRods{1}(2).Phi_h(Xs(i), obj.robot_linkage.CVRods{1}(2).Phi_dof, obj.robot_linkage.CVRods{1}(2).Phi_odr);
    
                    % Actuation Matrix
                    [Btau, ~] = distributedActuationMatrix(obj, q, Xs(i));
                    
                    % Actuator Path
                    A(obj.joint_dof + 1:obj.joint_dof + obj.strain_dof, :) = A + Ws(i)*(Bq')*Btau;
                end
            end

            % Compute QR Factorization
            if nargout > 1
                % A'P = QR
                [~, ~, P] = qr(A');

                % Compute Aa and Au
                if nargout > 2
                    Ap = (P')*A;
                    Aa = Ap(1:obj.m, :);
                    Au = Ap(obj.m + 1:end, :);
                end
            end
        end
        
        %% Test Actuation Matrix
        function testActuation(obj, options)
            arguments
                obj
                options.N = 1000;
            end
            % Candidates
            q_test = [zeros(obj.n, 1), randn(obj.n, options.N)];
            obj.q_sing = [];

            % Test the candidates
            for i = 1:(options.N + 1)
                % Store if A is singular
                if(rank(obj.actuationMatrix(q_test(:, i))) ~= obj.m)
                    obj.q_sing = [obj.q_sing, q_test(:, i)];
                end
            end

            % Warning
            if(~isempty(obj.q_sing))
                warning("The Actuation Matrix is singular for the following values:");
                disp(obj.q_sing);
            else
                disp("No singularities found in the Actuation Matrix.")
            end
        end

        %% Actuators Length
        function La = actuatorLengths(obj, q)
            % Init
            La = zeros(obj.n_sact, 1);
        
            % Compute only if there exists soft actuators
            if obj.n_sact > 0
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
        end

        function [theta, theta_dot] = transform(obj, q, qdot)
            % This function transform the coordinate q in the
            % collocated theta.

            % For now only underactuated system
            assert(obj.m <= obj.n, "Class supports only Underactuated system for now.");
            
            % To apply Pustina's Theorem, we need to reorder the rows of A
            % in such a way, the first m rows are full rank.
            % To perform this automatically, we use QR algorithm.

            % Compute Actuation Matrix
            [A, P] = obj.actuationMatrix(q);
            assert(rank(A) == obj.m, "Actuation Matrix is not full rank.");

            % Passive Output
            y = zeros(obj.m, 1);

            % Insert 1D Joint
            if obj.joint_dof > 0
                y(1) = q(1:obj.joint_dof);
            end

            if obj.n_sact > 0
                y(obj.joint_dof + 1: obj.joint_dof + obj.n_sact) = obj.actuatorLengths(q);
            end

            % Change of Coordinates for Underactuated system
            qp = (P')*q;            % qa = qp(1:obj.m) | qu = q(obj.m + 1:end)
            qdotp = (P')*qdot;      % qdota = qdotp(1:obj.m) | qdotu = qdot(obj.m + 1:end)

            % Formal Expression
            % theta = [y; zeros(obj.n - obj.m, 1)] + blkdiag(zeros(obj.m, obj.m), eye(obj.n - obj.m))*qp;
            % theta_dot = [(A')*qdot; [zeros(obj.n - obj.m, obj.m), eye(obj.n - obj.m)]*qdotp];

            % More efficient (and natural) expressions
            % Notice that should be (Ap')*qdotp, but = (A'*P)*P'*qdot = A'*qdot
            theta = [y; qp((obj.m + 1):end)];
            theta_dot = [(A')*qdot; qdotp((obj.m + 1):end)];
        end

        function z = groupTransform(obj, x)
            [theta, theta_dot] = obj.transform(x(1:obj.n), x(obj.n + 1:end));
            z = [theta; theta_dot];
        end

        function [q, qdot] = inverseTransformation(obj, theta, theta_dot, options)
            arguments
                obj
                theta
                theta_dot
                options.initial_guess = zeros(2*obj.n, 1);
                options.iter = 'iter';
            end
            % Numerical Solution
            tf_handle = @(x) [theta; theta_dot] - obj.groupTransform(x);

            % Set fsolve
            fsolve_opt = optimoptions('fsolve', 'Display', options.iter);
            x = fsolve(tf_handle, options.initial_guess, fsolve_opt);

            % Get Solution
            q = x(1:obj.n);
            qdot = x(obj.n + 1:end);
        end


        function [M, G, K, D] = dynamicMarices(obj, q, qdot)
            [~, dynamic_matrices] = obj.robot_linkage.my_dynamicsSolver(0, [q; qdot], zeros(obj.m, 1));
            % Inertia
            M = dynamic_matrices.M;
            G = - dynamic_matrices.G;
            % Elasticity vector
            K = dynamic_matrices.K*q;
            D = dynamic_matrices.D;
        end
    
        function [M_theta, G_theta, K_theta, D_theta] = transformSystem(obj, q, qdot)
            % This function transform the system in the collocated form.
            arguments
                obj
                q
                qdot
            end

            %% Compute Actuation Matrix
            [~, ~, Aa, Au] = obj.actuationMatrix(q);

            % Precomputing useful variables
            AaT = Aa';
            AuT = Au';
            invAaT = pinv(AaT);

            %% Compute Jacobian of the Transf. of Coordinates
            Jh = [AaT, AuT; zeros(obj.n - obj.m, obj.m), eye(obj.n - obj.m)];
            invJh = [invAaT, -invAaT*AuT; zeros(obj.n - obj.m, obj.m), eye(obj.n - obj.m)];

            %% Applying transformation to the Matrices (for now only G and K)
            [~, dynamic_matrices] = obj.robot_linkage.my_dynamicsSolver(0, [q; qdot], zeros(obj.m, 1));

            % Inertia
            M_theta = (invJh')*dynamic_matrices.M*invJh;
            G_theta = (invJh')*dynamic_matrices.G;
            % Elasticity vector
            K_theta = (invJh')*dynamic_matrices.K*q;
            D_theta = (invJh')*dynamic_matrices.D*invJh;
        end
    end
end