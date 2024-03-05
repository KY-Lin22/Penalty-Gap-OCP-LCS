function [z_Opt, Info] = solve_NLP_single(self, z_Init, p)
% solving only single NLP with given p
% NLP has the form:
%  min  J(z, p),
%  s.t. h(z, p) = 0,
% where: z is the variable,
%        p is the parameters,
%        J is the cost, and h is the constraint arranged in a stagewise manner.
%
% Syntax:
%          [z_Opt, Info] = solve_NLP_single(self, z_Init, p)
%          [z_Opt, Info] = self.solve_NLP_single(z_Init, p)
% Argument:
%          z_Init: double, NLP.Dim.z X 1, initial guess
%          p: double, given parameter
% Output:
%          z_Opt: double, NLP.Dim.z X 1, optimal solution found by solver
%          Info: struct, record the iteration information
%% check input
% check input z_Init and p
if ~all(size(z_Init) == [self.NLP.Dim.z, 1])
    error('z_Init has wrong dimension')
end
if ~all(size(p) == [self.NLP.Dim.p, 1])
    error('p has wrong dimension')
end
% check parameter
if (p(1) < 0)
    error('penalty parameter mu (i.e, p_1) should be nonnegative')
end

%% iteration routine (z: previous iterate z_{k-1}, z_k: current iterate z_{k}) 
% time record
Time = struct('gradEval', 0, 'searchDirection', 0, 'lineSearch', 0, 'else', 0, 'total', 0);
% load constant matrix
h_grad = self.NLP.FuncObj.h_grad; 
KKT_matrix_constant = self.KKT_matrix_constant;

k = 1;
beta = self.Option.LineSearch.betaInit;
z = z_Init;
gamma_h = zeros(self.NLP.Dim.h, 1);
J = full(self.NLP.FuncObj.J(z, p));
h = full(self.NLP.FuncObj.h(z, p));
while true
    %% step 0: check iteration counter
    if k > self.Option.maxIterNum
        % failure case 1: exceed the maximum number of iteration
        terminalStatus = 0;
        terminalMsg = ['- Solver fails: ', 'because the maximum number of iteration exceeded'];
        break
    end
    timeStart_total = tic; 

    %% step 1: Jacobian, Hessian, KKT error and gap evaluation of previous iterate z
    timeStart_gradEval = tic;

    % cost Jacobian
    J_grad = full(self.NLP.FuncObj.J_grad(z, p));
    % penalty hessian 
    switch self.Option.penalty_hessian_regularization
        case 0
            % without regularization
            J_penalty_hessian = sparse(self.NLP.FuncObj.J_penalty_hessian(z, p));
        case 1
            % with regularization
            w = full(self.NLP.FuncObj.w(z));
            J_penalty_hessian = sparse(self.NLP.FuncObj.J_penalty_hessian_regular(z, p, w)) ...
                + self.Option.penalty_hessian_additional_regular_param * speye(self.NLP.Dim.z);
    end
    % KKT error (L_inf norm)
    LAG_grad_z = J_grad + gamma_h' * h_grad;
    KKT_error_primal = norm(h, inf);
    KKT_error_dual = norm(LAG_grad_z, inf);

    timeElasped_gradEval = toc(timeStart_gradEval);

    %% step 2: search direction evaluation based on previous iterate z
    timeStart_searchDirection = tic;

    [dz, gamma_h_k, Info_SearchDirection] = self.evaluate_search_direction(h, J_grad, J_penalty_hessian);

    % check status
    if Info_SearchDirection.status == 0
        % failure case 2: qp solver fails
        terminalStatus = -1;
        terminalMsg = ['- Solver fails: ', Info_SearchDirection.msg];
        break
    else
        % dzNorm (L_inf norm)
        dzNorm = norm(dz, inf);
    end

    timeElasped_searchDirection = toc(timeStart_searchDirection);

    %% step 3: check whether we can terminate successfully based on the previous iterate z
    if max([KKT_error_primal, KKT_error_dual]) < self.Option.tol.KKT_error_total
        % Success case 1: the KKT error satisfies tolerance
        terminalStatus = 1;
        terminalMsg = ['- Solver succeeds: ', 'because the KKT error satisfies tolerance'];
        break
    elseif dzNorm < self.Option.tol.dzNorm
        % Success case 2: the norm of search direction satisfies tolerance
        terminalStatus = 1;
        terminalMsg = ['- Solver succeeds: ', 'because the norm of search direction satisfies tolerance'];      
        break
    elseif (KKT_error_primal <= self.Option.tol.KKT_error_primal) && (KKT_error_dual <= self.Option.tol.KKT_error_dual)
        % Success case 3: primal and dual error satisfy tolerance
        terminalStatus = 1;
        terminalMsg = ['- Solver succeeds: ', 'because primal and dual error satisfy tolerance']; 
        break
    end  

    %% step 4: merit line search
    timeStart_lineSearch = tic;

    [z_k, Info_LineSearch] = self.line_search_merit(beta, z, dz, p, J, h, J_grad, J_penalty_hessian);
    % check status
    if Info_LineSearch.status == 0
        % failure case 3: line search fails
        terminalStatus = -2;
        terminalMsg = ['- Solver fails: ', 'because merit line search reaches the min stepsize'];        
        break
    else
        % extract quantities (J, h) associated with z_k
        J_k      = Info_LineSearch.J;
        h_k      = Info_LineSearch.h;
        beta_k   = Info_LineSearch.beta;
        stepSize = Info_LineSearch.stepSize;
        merit    = Info_LineSearch.merit;
    end  
    timeElasped_lineSearch = toc(timeStart_lineSearch);

    %% step 5: record and print information of this iteration k
    timeElasped_total = toc(timeStart_total);

    Time.gradEval = Time.gradEval + timeElasped_gradEval;
    Time.searchDirection = Time.searchDirection + timeElasped_searchDirection; 
    Time.lineSearch = Time.lineSearch + timeElasped_lineSearch;
    Time.total = Time.total + timeElasped_total;
    % print
    if self.Option.printLevel == 2
        % head
        if mod(k, 10) == 1
            disp('----------------------------------------------------------------------------------------')
            headMsg = ' Iter | cost(ocp/penalty) |  gapRes  | KKT(primal/dual)|  dzNorm  |   beta   | stepsize |  merit   | merit(t) | time(ms) |';
            disp(headMsg)
        end
        % previous iterate message
        prevIterMsg = [' ',...
            num2str(k,'%10.3d'),'  | ',...
            num2str(full(self.NLP.FuncObj.J_ocp(z, p)), '%10.2e'), ' ',...
            num2str(full(self.NLP.FuncObj.J_penalty(z, p)), '%10.2e'),' | ',...
            num2str(norm(full(self.NLP.FuncObj.D_gap_func(z)), inf), '%10.2e'), ' | ',...
            num2str(KKT_error_primal, '%10.1e'), ' ',...
            num2str(KKT_error_dual, '%10.1e'),' | ',...
            num2str(dzNorm,'%10.2e'), ' | ',...
            num2str(beta_k,'%10.2e'), ' | ',...
            num2str(stepSize,'%10.2e'), ' | ',...
            num2str(merit(1),'%10.2e'), ' | ',...
            num2str(merit(2),'%10.2e'), ' | ',...
            num2str(1000 * timeElasped_total,'%10.2e'), ' | '];
        disp(prevIterMsg)
    end

    %% step 6: prepare next iteration
    k = k + 1;
    beta = beta_k;
    z = z_k;
    gamma_h = gamma_h + stepSize * (gamma_h_k - gamma_h);
    J = J_k;
    h = h_k;
end

%% return optimal solution and create information
% return previous iterate as solution
z_Opt = z;
% create Info (basic: time, iterNum, terminalStatus)
Time.else = Time.total - Time.gradEval - Time.searchDirection - Time.lineSearch;
Info.Time = Time;
Info.iterNum = k - 1;
Info.terminalStatus = terminalStatus;
Info.terminalMsg = terminalMsg;
% create Info (corresponds to the solution z_Opt: dual variable, cost, KKT, natural residual)
Info.gamma_h             = gamma_h;
Info.cost_ocp            = full(self.NLP.FuncObj.J_ocp(z, p));
Info.cost_penalty        = full(self.NLP.FuncObj.J_penalty(z, p));
Info.KKT_error_primal    = KKT_error_primal;
Info.KKT_error_dual      = KKT_error_dual;
Info.VI_natural_residual = self.evaluate_natural_residual(z);
% display termination and solution message
if (self.Option.printLevel == 1) || (self.Option.printLevel == 2)
    disp('*------------------ Solution Information ------------------*')
    disp('1. Terminal Status')
    disp(Info.terminalMsg)
    disp('2. Iteration Process Message')
    disp(['- Iterations: ................', num2str(Info.iterNum)])
    disp(['- TimeElapsed: ...............', num2str(Info.Time.total,'%10.3f'), 's'])
    disp(['- AverageTime: ...............', num2str(1000 * Info.Time.total /Info.iterNum,'%10.2f'), ' ms/Iter'])
    disp('3. Solution Message')
    disp(['- Cost (ocp): ................', num2str(Info.cost_ocp,'%10.3e'), '; '])
    disp(['- Cost (penalty): ............', num2str(Info.cost_penalty,'%10.3e'), '; '])
    disp(['- KKT (primal): ..............', num2str(Info.KKT_error_primal,'%10.3e'), '; '])
    disp(['- KKT (dual): ................', num2str(Info.KKT_error_dual,'%10.3e'), '; '])
    disp(['- VI natural residual: .......', num2str(Info.VI_natural_residual,'%10.3e'), '; '])
end

end