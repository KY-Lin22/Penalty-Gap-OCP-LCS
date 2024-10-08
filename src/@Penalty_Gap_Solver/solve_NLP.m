function [z_Opt, Info] = solve_NLP(self, z_Init, p_Init, p_End)
% solve NLP with given z_Init, p_Init, and p_End by gap_penalty solver using continuation method
% NLP has the form:
%  min  J(z, p),
%  s.t. h(z, p) = 0,
% where: z is the variable,
%        p is the parameters,
%        J is the cost, and h is the constraints
% Syntax:
%          [z_Opt, Info] = solve_NLP(self, z_Init, p_Init, p_End)
%          [z_Opt, Info] = self.solve_NLP(z_Init, p_Init, p_End)
% Argument:
%          z_Init: double, NLP.Dim.z X 1, initial guess
%          p_Init: double, problem parameter (initial) p = mu
%          p_End: double, problem parameter (end) p = mu
% Output:
%          z_Opt: double, NLP.Dim.z X 1, optimal solution found by solver
%          Info: struct, record the iteration information

import casadi.*

%% check option and input
% check option (TODO)
self.Option.printLevel = 0; % do not print each homotopy problem iteration information
% check input z_Init
if ~all(size(z_Init) == [self.NLP.Dim.z, 1])
    error('z_Init has wrong dimension')
end
% check parameter
if ~all(size(p_Init) == [self.NLP.Dim.p, 1])
    error('p_Init has wrong dimension')
end
if ~all(size(p_End) == [self.NLP.Dim.p, 1])
    error('p_End has wrong dimension')
end
% check penalty parameter
if (p_Init(1) < 0) || (p_End(1) < 0)
    error('penalty parameter mu (i.e., p_1) should be nonnegative')
end
if p_Init(1) > p_End(1)
    error('mu_Init should not larger than mu_End')
end

%% create record for time and log 
% evaluate the max number of continuation step based on given mu_Init, mu_End 
% and specified kappa_mu_times
mu_Init = p_Init(1);
mu_End = p_End(1);
mu_test = mu_Init;
continuationStepMaxNum = 1;
while true
    if mu_test == mu_End
        break
    else        
        mu_test = min([self.Option.Homotopy.kappa_mu_times * mu_test, mu_End]);
        continuationStepMaxNum = continuationStepMaxNum + 1;
    end
end
% record
Time = struct('gradEval', 0, 'searchDirection', 0, 'lineSearch', 0, 'else', 0, 'total', 0);
iterNum = 0;

%% continuation loop (j: continuation step counter)
% polish z_Init
z_Init = self.polish_initial_guess(z_Init);
z_Init_j = z_Init;
p_j = p_Init; 
j = 1;
while true
    %% step 1: solve a NLP with given p
    [z_Opt_j, Info_j] = self.solve_NLP_single(z_Init_j, p_j);
    dual_var_Opt_j     = Info_j.gamma_h;
    J_ocp_j            = Info_j.cost_ocp;
    J_penalty_j        = Info_j.cost_penalty;
    KKT_error_primal_j = Info_j.KKT_error_primal;
    KKT_error_dual_j   = Info_j.KKT_error_dual; 
    VI_nat_res_j       = Info_j.VI_natural_residual;
    iterNum_j          = Info_j.iterNum;

    %% step 2: record and print information of the current continuation iterate
    Time.gradEval        = Time.gradEval        + Info_j.Time.gradEval;
    Time.searchDirection = Time.searchDirection + Info_j.Time.searchDirection;
    Time.lineSearch      = Time.lineSearch      + Info_j.Time.lineSearch;
    Time.else            = Time.else            + Info_j.Time.else;
    Time.total           = Time.total           + Info_j.Time.total;  
    iterNum              = iterNum              + iterNum_j;
    if mod(j, 10) == 1
        disp('---------------------------------------------------------------------------------------------')
        headMsg = ' step  | param(mu)| cost(ocp/penalty) | KKT(primal/dual)| nat_res | iterNum | time(s) ';
        disp(headMsg)
    end
    prevIterMsg = [' ',...
        num2str(j,'%10.2d'), '/', num2str(continuationStepMaxNum,'%10.2d'),' |  ',...
        num2str(p_j(1), '%10.1e'), ' | ',...
        num2str(J_ocp_j, '%10.2e'), ' ', num2str(J_penalty_j, '%10.2e'),' | ',...
        num2str(KKT_error_primal_j, '%10.1e'), ' ', num2str(KKT_error_dual_j, '%10.1e'),' | ',...
        num2str(VI_nat_res_j, '%10.1e'),' |   ',...
        num2str(iterNum_j, '%10.4d'),'  | ',...
        num2str(Info_j.Time.total, '%10.4f')];
    disp(prevIterMsg)

    %% step 3: check ternimation based on the current homotopy iterate
    if (Info_j.terminalStatus == 1) && (VI_nat_res_j <= self.Option.Homotopy.VI_nat_res_tol)
        % IPOPT at this homotopy iteration finds the optimal solution satisfying the desired VI natural residual
        terminalStatus = 1;
        terminalMsg = Info_j.terminalMsg;
        break
    elseif (Info_j.terminalStatus ~= 1)
        % IPOPT at this homotopy iteration fails to find the optimal solution
        terminalStatus = 0;
        terminalMsg = Info_j.terminalMsg;
        break
    elseif j == continuationStepMaxNum
        % IPOPT still can not find the optimal solution in the final homotopy iteration
        terminalStatus = 0;
        terminalMsg = 'solver can not find the optimal solution satisfying the desired VI natural residual';
        break
    else
        % IPOPT at this homotopy iteration (not the final) finds the optimal solution, prepare for next homotopy
        % update initial guess
        z_Init_j = z_Opt_j;
        % update penalty parameter
        mu_j = p_j(1);
        mu_j = min([self.Option.Homotopy.kappa_mu_times * mu_j, mu_End]);
        % update parameter vector
        p_j(1) = mu_j;
        % update continuation step counter
        j = j + 1;
    end

end

%% reshape and return the current homotopy iterate as the optimal solution
Z_Opt_j = reshape(z_Opt_j, self.NLP.Dim.z_Node(end), self.OCP.nStages);

X_Opt_j  = Z_Opt_j(1 : self.NLP.Dim.z_Node(1), :);
U_Opt_j  = Z_Opt_j(self.NLP.Dim.z_Node(1) + 1 : self.NLP.Dim.z_Node(2), :);
XI_Opt_j = Z_Opt_j(self.NLP.Dim.z_Node(2) + 1 : self.NLP.Dim.z_Node(4), :);

[~, decouple_func_lambda, decouple_func_eta] = ...
    self.NLP.create_element_wise_concatenation_func(self.OCP);

decouple_func_lambda_map = decouple_func_lambda.map(self.OCP.nStages);
decouple_func_eta_map = decouple_func_eta.map(self.OCP.nStages);
LAMBDA_Opt_j = full(decouple_func_lambda_map(XI_Opt_j));
ETA_Opt_j = full(decouple_func_eta_map(XI_Opt_j));
z_Opt_j = reshape([X_Opt_j; U_Opt_j; LAMBDA_Opt_j; ETA_Opt_j],...
    (self.OCP.Dim.x + self.OCP.Dim.u + 2*self.OCP.Dim.lambda) * self.OCP.nStages, 1);

z_Opt = z_Opt_j;
% create Info
Info.continuationStepNum = j;
Info.terminalStatus = terminalStatus;
Info.terminalMsg = terminalMsg;
Info.dual_var = dual_var_Opt_j;
Info.cost.ocp = J_ocp_j;
Info.cost.penalty = J_penalty_j;
Info.KKT_error.primal = KKT_error_primal_j;
Info.KKT_error.dual = KKT_error_dual_j;
Info.VI_natural_residual = VI_nat_res_j;
Info.Time = Time;
Info.time = Time.total;
Info.iterNum = iterNum;
% display homotopy terminal result and then break
disp('*----------------------------------- Solution Information ------------------------------------*')
disp(['1. Terminal Message: ', Info.terminalMsg])
disp('2. Continuation Step Message')
disp(['- TimeElapsed: ................... ', num2str(Info.time,'%10.3f'), ' s'])
disp(['- Continuation Step: ............. ', num2str(Info.continuationStepNum)])
disp(['- Time Per Continuation Step: .... ', num2str(Info.time / Info.continuationStepNum,'%10.3f'), ' s/Step'])
disp(['- Iterations: .................... ', num2str(Info.iterNum)])
disp(['- Time Per Iteration: ............ ', num2str(1000 * Info.time / Info.iterNum,'%10.3f'), ' ms/Iter'])
disp('3. Solution Message')
disp(['- Cost(ocp): ..................... ', num2str(Info.cost.ocp,'%10.3e'), '; '])
disp(['- Cost(penalty): ................. ', num2str(Info.cost.penalty,'%10.3e'), '; '])
disp(['- KKT(primal): ................... ', num2str(Info.KKT_error.primal,'%10.3e'), '; '])
disp(['- KKT(dual): ..................... ', num2str(Info.KKT_error.dual,'%10.3e')  '; '])
disp(['- natural residual: .............. ', num2str(Info.VI_natural_residual,'%10.3e'), '; '])
end