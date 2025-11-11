function Gamma = potential_hessian(robot_linkage, q, qd, u)
    %% Find q2dot
    t = 0; % Time-invariant

    % dynamicsSolver function
    [y,~,~,~] = robot_linkage.dynamicsSolver(t, [q; qd], u);
    
    % dxdt = [q_dot; q_2dot]
    qdd =  y(1:robot_linkage.ndof);

    %% DAE Jacobian
    dae_index = 1;
    % lambda is required in case of closed-loop joint
    lambda = [];

    % DAE Jacobian
    [~,~,~,dID_dq,~,~,~,dtau_dq,~,~,~,~,~,~,~] = robot_linkage.DAEJacobians(t, q, qd, qdd, u, lambda, dae_index);

    % Compute the Hessian
    Gamma = dtau_dq - dID_dq;
end