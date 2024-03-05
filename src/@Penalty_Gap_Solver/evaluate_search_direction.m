function [dz, gamma_h_k, Info] = evaluate_search_direction(self, h, J_grad, J_penalty_hessian)
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

%% solve qp
switch self.Option.qpSolver
    case 'fbstab_sparse'
        % set QP problem (follow the NLP KKT stationary condition) and solver option
        LAG_hessian = self.NLP.FuncObj.J_ocp_hessian + J_penalty_hessian;
        qp_Prob = struct('H', LAG_hessian, 'f', J_grad',...
            'G', self.NLP.FuncObj.h_grad, 'h', -h, 'A', sparse(1, self.NLP.Dim.z), 'b', ones(1, 1));
        % solve QP and extract solution
        [x, out] = fbstab_sparse(self.Option.fbstab_x0, qp_Prob, self.Option.fbstab_options);
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
        LAG_hessian = self.NLP.FuncObj.J_ocp_hessian + J_penalty_hessian;
        % set QP problem (follow the NLP KKT stationary condition) and solver option   
        self.QP_solver = osqp;
        self.QP_solver.setup(LAG_hessian, J_grad',...
            self.NLP.FuncObj.h_grad, -h, -h, self.Option.osqp_options)% P, q, A, l, u
        % initial guess
        self.QP_solver.warm_start('x', self.Option.osqp_x0, 'y', self.Option.osqp_y0);
        % solve QP and extract solution
        results = self.QP_solver.solve();
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
        % KKT residual
        KKT_residual = [-J_grad'; -h];
        % KKT matrix
        [i_J_pen_hess, j_J_pen_hess, s_J_pen_hess] = find(J_penalty_hessian);
        KKT_matrix_update = sparse(i_J_pen_hess, j_J_pen_hess, s_J_pen_hess,...
            self.NLP.Dim.z + self.NLP.Dim.h, self.NLP.Dim.z + self.NLP.Dim.h, length(s_J_pen_hess));
        KKT_matrix = self.KKT_matrix_constant + KKT_matrix_update;
        % solve linear system
        dz_gamma_h_k = KKT_matrix\KKT_residual;
        dz = dz_gamma_h_k(1 : self.NLP.Dim.z, 1);
        gamma_h_k = dz_gamma_h_k(self.NLP.Dim.z + 1 : end, 1);
        QPstatus = 1;
        QPmsg = 'succeeds';
end

%% create info
Info.status = QPstatus;
Info.msg = QPmsg;

end