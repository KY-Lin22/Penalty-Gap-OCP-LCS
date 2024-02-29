function [dz, gamma_h_k, Info] = evaluate_search_direction(self, h, J_grad, h_grad, LAG_hessian)
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
        % initial guess
        x0.z = zeros(NLP.Dim.z, 1); % primal variable
        x0.l = zeros(NLP.Dim.h, 1); % dual variable for equality constraint h
        x0.v = zeros(1, 1); % dual variable for inequality constraint c
        % set QP problem (follow the NLP KKT stationary condition) and solver option
        qp_Prob = struct('H', LAG_hessian, 'f', J_grad',...
            'G', h_grad, 'h', -h, 'A', sparse(1, NLP.Dim.z), 'b', ones(1, 1));
        qp_Option = Option.fbstab_options;

        % solve QP and extract solution
        timeStart = tic;
        [x, out] = fbstab_sparse(x0, qp_Prob, qp_Option);
        timeElasped = toc(timeStart);
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
        qp_Prob = struct('P', LAG_hessian, 'q', J_grad',...
            'A', h_grad, 'l', -h, 'u', -h);       
        qp_Option = Option.osqp_options;
        osqpSolver.setup(qp_Prob.P, qp_Prob.q, qp_Prob.A, qp_Prob.l, qp_Prob.u, qp_Option)

        % initial guess
        x0 = zeros(NLP.Dim.z, 1); % primal variable
        y0 = zeros(NLP.Dim.h, 1); % dual variable for equality constraint h and inequality constraint c 
        osqpSolver.warm_start('x', x0, 'y', y0);

        % solve QP and extract solution
        timeStart = tic;
        results = osqpSolver.solve();
        timeElasped = toc(timeStart);
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
end

%% create info
Info.time = timeElasped;
Info.status = QPstatus;
Info.msg = QPmsg;

end