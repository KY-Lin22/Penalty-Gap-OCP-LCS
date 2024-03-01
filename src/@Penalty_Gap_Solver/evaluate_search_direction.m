function [dz, gamma_h_k, Info] = evaluate_search_direction(self, h, J_grad, h_grad, J_ocp_hessian, J_penalty_hessian)
%evaluate the search direction by solving a QP problem
%
% fbstab_sparse solver requires the QP problem has the form:
% min.  1/2 z'Hz + f'z
% s.t.  Gz = h
%       Az <= b
%
% osqp solver requires the QP problem has the form:
% min.  1/2 x'Px + q'x
% s.t.  l <= Ax <= u

NLP = self.NLP;
Option = self.Option;

%% solve qp
switch Option.qpSolver
    case 'fbstab_sparse'
        % set QP problem (follow the NLP KKT stationary condition) and solver option
        qp_Prob = struct('H', J_ocp_hessian + J_penalty_hessian, 'f', J_grad',...
            'G', h_grad, 'h', -h, 'A', sparse(1, NLP.Dim.z), 'b', ones(1, 1));
        % solve QP and extract solution
        [x, out] = fbstab_sparse(Option.fbstab_x0, qp_Prob, Option.fbstab_options);
        dz = x.z; % optimal primal variable
        gamma_h_k = x.l; % optimal dual variable for equality constraint h
        % return solver status
        if out.eflag == 0
            % QP solver finds the optimal solution
            QPstatus = 1;
            QPmsg = 'fbstab_sparse succeeds';
        else
            % QP solver fails to find the optimal solution
            QPstatus = 0;
            switch out.eflag
                case -1
                    QPmsg = ['fbstab_sparse fails because ', 'Maximum number of iterations exceeded'];
                case -2
                    QPmsg = ['fbstab_sparse fails because ', 'Algorithm stalled'];
                case -3
                    QPmsg = ['fbstab_sparse fails because ', 'Problem is infeasible'];
                case -4
                    QPmsg = ['fbstab_sparse fails because ', 'Problem is unbounded below (dual infeasible)'];
                case -5
                    QPmsg = ['fbstab_sparse fails because ', 'Problem is primal and dual infeasible'];
            end
        end

    case 'osqp'
        % initialization
        osqpSolver = osqp;
        % set QP problem (follow the NLP KKT stationary condition) and solver option
        qp_Prob = struct('P', J_ocp_hessian + J_penalty_hessian, 'q', J_grad',...
            'A', h_grad, 'l', -h, 'u', -h);       
        osqpSolver.setup(qp_Prob.P, qp_Prob.q, qp_Prob.A, qp_Prob.l, qp_Prob.u, Option.osqp_options)
        % initial guess
        osqpSolver.warm_start('x', Option.osqp_x0, 'y', Option.osqp_y0);
        % solve QP and extract solution
        results = osqpSolver.solve();
        dz = results.x; % optimal primal variable
        gamma_h_k = results.y; % optimal dual variable for equality constraint h
        % return solver status
        if results.info.status_val == 1
            % QP solver finds the optimal solution
            QPstatus = 1;
            QPmsg = 'osqp succeeds';
        else
            % QP solver fails to find the optimal solution
            QPstatus = 0;
            QPmsg = ['osqp fails because ', results.info.status] ;
        end

    case 'direct_sparse'
        % only suitable when does not have inequality constraint
        KKT_matrix_constant = self.KKT_matrix_constant;
        % KKT residual
        KKT_residual = [-J_grad'; -h];
        % KKT matrix
        [i_J_pen_hess, j_J_pen_hess, s_J_pen_hess] = find(J_penalty_hessian);
        KKT_matrix_update = sparse(i_J_pen_hess, j_J_pen_hess, s_J_pen_hess,...
            NLP.Dim.z + NLP.Dim.h, NLP.Dim.z + NLP.Dim.h, length(s_J_pen_hess));
        KKT_matrix = KKT_matrix_constant + KKT_matrix_update;
        % solve linear system
        dz_gamma_h_k = KKT_matrix\KKT_residual;
        dz = dz_gamma_h_k(1 : NLP.Dim.z, 1);
        gamma_h_k = dz_gamma_h_k(NLP.Dim.z + 1 : end, 1);
        QPstatus = 1;
        QPmsg = 'succeeds';
end

%% create info
Info.status = QPstatus;
Info.msg = QPmsg;

end