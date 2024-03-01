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
%% check option and input
% check option (TODO: write a function)
if self.Option.printLevel == 2
    self.Option.recordLevel = 1;
end
NLP = self.NLP;
Option = self.Option;
% check input z_Init
if ~all(size(z_Init) == [NLP.Dim.z, 1])
    error('z_Init has wrong dimension')
end
% check parameter
if (p(1) < 0)
    error('penalty parameter mu (i.e, p_1) should be nonnegative')
end

%% create record for time and log
% time
Time = struct('gradEval', 0, 'searchDirection', 0, 'lineSearch', 0, 'else', 0, 'total', 0);
% log (also record some quantities that may not be used in the iteration routine,
% e.g., J_ocp, J_penalty, D_gap_func_res, D_gap_grad_res)
if Option.recordLevel == 1
    Log.cost      = zeros(Option.maxIterNum, 3); % [ocp, penalty, total]
    Log.gap       = zeros(Option.maxIterNum, 2); % [max D_gap_func_res, max D_gap_grad_res]
    Log.KKT_error = zeros(Option.maxIterNum, 3); % [primal, dual, total]   
    Log.dzNorm    = zeros(Option.maxIterNum, 1);
    Log.beta      = zeros(Option.maxIterNum, 1); 
    Log.stepSize  = zeros(Option.maxIterNum, 1);
    Log.merit     = zeros(Option.maxIterNum, 2); % [merit, merit_k]
end

%% prepare the first iteration (z: previous iterate z_{k-1}, z_k: current iterate z_{k}) 
if (Option.printLevel == 1) || (Option.printLevel == 2)
    disp('**************************************************************************************');
    disp('Initializing...');
end

% counter, beta, z, multipler, cost and constraint function
k = 1;
beta = Option.LineSearch.betaInit;
z = z_Init;
gamma_h = zeros(NLP.Dim.h, 1);
J = full(NLP.FuncObj.J(z, p));
h = full(NLP.FuncObj.h(z, p));

% constant matrix
h_grad = NLP.FuncObj.h_grad;
J_ocp_hessian = NLP.FuncObj.J_ocp_hessian;

%% iteration routine 
% begin iteration routine
if (Option.printLevel == 1) || (Option.printLevel == 2)
    disp('Computing the Optimal Solution...');
end
% iteration loop
while true
    %% step 0: check iteration counter
    if k > Option.maxIterNum
        % failure case 1: exceed the maximum number of iteration
        terminalStatus = 0;
        terminalMsg = ['- Solver fails: ', 'because the maximum number of iteration exceeded'];
        break
    end
    timeStart_total = tic; 

    %% step 1: Jacobian, Hessian, KKT error and gap evaluation of previous iterate z
    timeStart_gradEval = tic;

    % cost Jacobian
    J_grad = full(NLP.FuncObj.J_grad(z, p));
    % D gap function and Jacobian
    D_gap_grad = full(NLP.FuncObj.D_gap_grad(z));
    D_gap_hessian = sparse(NLP.FuncObj.D_gap_hessian(z));
    % Huber hessian
    Huber_hessian = sparse(NLP.Huber_hessian(D_gap_grad));
    % penalty Hessian
    J_penalty_hessian = D_gap_hessian' * Huber_hessian * D_gap_hessian;

    % KKT error (L_inf norm)
    KKT_error_primal = norm(h, inf);
    KKT_error_dual = norm(J_grad' + h_grad' * gamma_h, inf);
    KKT_error_total = max([KKT_error_primal, KKT_error_dual]);

    timeElasped_gradEval = toc(timeStart_gradEval);

    %% step 2: search direction evaluation based on previous iterate z
    timeStart_searchDirection = tic;

    % solving a sparse QP subproblem    
    [dz, gamma_h_k, Info_SearchDirection] = ...
        self.evaluate_search_direction(h, J_grad, h_grad, J_ocp_hessian, J_penalty_hessian);
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
    if KKT_error_total < Option.tol.KKT_error_total
        % Success case 1: the KKT error satisfies tolerance
        terminalStatus = 1;
        terminalMsg = ['- Solver succeeds: ', 'because the KKT error satisfies tolerance'];
        break
    elseif dzNorm < Option.tol.dzNorm
        % Success case 2: the norm of search direction satisfies tolerance
        terminalStatus = 1;
        terminalMsg = ['- Solver succeeds: ', 'because the norm of search direction satisfies tolerance'];      
        break
    elseif (KKT_error_primal <= Option.tol.KKT_error_primal) && (KKT_error_dual <= Option.tol.KKT_error_dual)
        % Success case 3: primal and dual error satisfy tolerance
        terminalStatus = 1;
        terminalMsg = ['- Solver succeeds: ', 'because primal and dual error satisfy tolerance']; 
        break
    end  

    %% step 4: merit line search
    timeStart_lineSearch = tic;

    [z_k, Info_LineSearch] = self.line_search_merit(beta, z, dz, p, J, h, J_grad, J_ocp_hessian, J_penalty_hessian);
    % check status
    if Info_LineSearch.status == 0
        % failure case 3: line search fails
        terminalStatus = -2;
        terminalMsg = ['- Solver fails: ', 'because merit line search reaches the min stepsize'];        
        break
    else
        % extract quantities (J, h) associated with z_k
        J_k = Info_LineSearch.J;
        h_k = Info_LineSearch.h;
        beta_k = Info_LineSearch.beta;
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

    if Option.recordLevel == 1
        % record (including some quantities that may not be used in the SGFL iteration rountie, e.g., J_ocp, J_penalty)
        Log.cost(k, :)      = [full(NLP.FuncObj.J_ocp(z, p)), full(NLP.FuncObj.J_penalty(z, p)), J];
        Log.gap(k, :)       = [norm(full(NLP.FuncObj.D_gap_func(z)), inf), norm(D_gap_grad, inf)];
        Log.KKT_error(k, :) = [KKT_error_primal, KKT_error_dual, KKT_error_total];        
        Log.dzNorm(k)       = dzNorm;
        Log.beta(k)         = beta_k;
        Log.stepSize(k)     = stepSize;    
        Log.merit(k, :)     = merit;
        % print
        if Option.printLevel == 2
            % head
            if mod(k, 10) == 1
                disp('----------------------------------------------------------------------------------------------------------------------------------------------------')
                headMsg = ' Iter | cost(ocp)| cost(pen)|  gapRes  |gapGradRes|  KKT(P)  |  KKT(D)  |  dzNorm  |   beta   | stepsize |  merit   | merit(t) | time(ms) |';
                disp(headMsg)
            end
            % previous iterate message
            prevIterMsg = [' ',...
                num2str(k,'%10.3d'),'  | ',...
                num2str(Log.cost(k, 1),'%10.2e'), ' | ',...
                num2str(Log.cost(k, 2),'%10.2e'), ' | ',...
                num2str(Log.gap(k, 1), '%10.2e'), ' | ',...
                num2str(Log.gap(k, 2), '%10.2e'), ' | ',...
                num2str(Log.KKT_error(k, 1), '%10.2e'), ' | ',...
                num2str(Log.KKT_error(k, 2), '%10.2e'), ' | ',...
                num2str(Log.dzNorm(k),'%10.2e'), ' | ',...
                num2str(Log.beta(k),'%10.2e'), ' | ',...
                num2str(Log.stepSize(k),'%10.2e'), ' | ',...
                num2str(Log.merit(k, 1),'%10.2e'), ' | ',...
                num2str(Log.merit(k, 2),'%10.2e'), ' | ',...
                num2str(1000 * timeElasped_total,'%10.2e'), ' | '];
            disp(prevIterMsg)
        end

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
Time.else = Time.total - Time.searchDirection - Time.lineSearch - Time.gradEval;
Info.Time = Time;
Info.iterNum = k - 1;
Info.terminalStatus = terminalStatus;
Info.terminalMsg = terminalMsg;
% create Info (corresponds to the solution z_Opt: dual variable, cost, KKT, natural residual)
Info.gamma_h                 = gamma_h;
Info.cost_ocp                = full(NLP.FuncObj.J_ocp(z, p));
Info.cost_penalty            = full(NLP.FuncObj.J_penalty(z, p));
Info.cost_total              = J;
Info.KKT_error_primal        = KKT_error_primal;
Info.KKT_error_dual          = KKT_error_dual;
Info.KKT_error_total         = KKT_error_total;
Info.VI_natural_residual     = self.evaluate_natural_residual(z);
% create Info (log)
if Option.recordLevel == 1
    Info.Log.cost      = Log.cost(1 : k - 1, :);
    Info.Log.gap       = Log.gap(1 : k - 1, :);
    Info.Log.KKT_error = Log.KKT_error(1 : k - 1, :);
    Info.Log.dzNorm    = Log.dzNorm(1 : k - 1, :);
    Info.Log.beta      = Log.beta(1 : k - 1, :);
    Info.Log.stepSize  = Log.stepSize(1 : k - 1, :);
    Info.Log.merit     = Log.merit(1 : k - 1, :);
end
% display termination and solution message, then break rountie
if (Option.printLevel == 1) || (Option.printLevel == 2)
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
    disp(['- Cost (total): ..............', num2str(Info.cost_total,'%10.3e'), '; '])
    disp(['- KKT (primal): ..............', num2str(Info.KKT_error_primal,'%10.3e'), '; '])
    disp(['- KKT (dual): ................', num2str(Info.KKT_error_dual,'%10.3e'), '; '])
    disp(['- KKT (total): ...............', num2str(Info.KKT_error_total,'%10.3e'), '; '])
    disp(['- VI natural residual: .......', num2str(Info.VI_natural_residual,'%10.3e'), '; '])
end

% end of the iteration routine 
if (Option.printLevel == 1) || (Option.printLevel == 2)
    disp('Done!')
    disp('**************************************************************************************');
end

end