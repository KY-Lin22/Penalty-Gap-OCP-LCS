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

NLP = self.NLP;
Option = self.Option;

% check input z_Init
if ~all(size(z_Init) == [NLP.Dim.z, 1])
    error('z_Init has wrong dimension')
end

% check penalty parameter
if (p_Init(1) < 0) || (p_End(1) < 0)
    error('penalty parameter mu (i.e., p_1) should be nonnegative')
end
if p_Init(1) > p_End(1)
    error('mu_Init should not larger than mu_End')
end
% load parameter
kappa_mu_times = Option.Homotopy.kappa_mu_times;
VI_nat_res_tol = Option.Homotopy.VI_nat_res_tol;

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
        mu_test = min([kappa_mu_times * mu_test, mu_End]);
        continuationStepMaxNum = continuationStepMaxNum + 1;
    end
end
% log (time, param, cost, KKT error, stepSize, natRes)
Time = struct('searchDirection', 0, 'lineSearch', 0, 'gradEval', 0, 'else', 0, 'total', 0);% detail of continuation step time
Log.param               = zeros(continuationStepMaxNum, 1); % mu 
Log.cost                = zeros(continuationStepMaxNum, 2); % [ocp, penalty]
Log.KKT_error           = zeros(continuationStepMaxNum, 2); % [primal, dual]
if Option.recordLevel == 1
    Log.stepSize_primal     = zeros(continuationStepMaxNum, 2); % [min, average]
    Log.stepSize_dual       = zeros(continuationStepMaxNum, 2); % [min, average]
end
Log.VI_natural_residual = zeros(continuationStepMaxNum, 1);
Log.iterNum             = zeros(continuationStepMaxNum, 1);
Log.timeElapsed         = zeros(continuationStepMaxNum, 1); % elapsed time in each continuation step

%% continuation loop (j: continuation step counter)
z_Init_j = z_Init;
p_j = p_Init; 
j = 1;
while true
    %% step 1: solve a NLP with given p
    [z_Opt_j, Info_j] = self.solve_NLP_single(z_Init_j, p_j);
    dual_var_Opt_j = full(Info_j.gamma_h);
    J_ocp_j = full(NLP.FuncObj.J_ocp(z_Opt_j, p_j));
    J_penalty_j = full(NLP.FuncObj.J_penalty(z_Opt_j, p_j));
    KKT_error_primal_j = Info_j.KKT_error_primal;
    KKT_error_dual_j = Info_j.KKT_error_dual; 
    VI_nat_res_j = Info_j.VI_natural_residual;

    %% step 2: record and print information of the current continuation iterate
    % time
    Time.searchDirection = Time.searchDirection + Info_j.Time.searchDirection;
    Time.lineSearch      = Time.lineSearch      + Info_j.Time.lineSearch;
    Time.gradEval        = Time.gradEval        + Info_j.Time.gradEval;
    Time.else            = Time.else            + Info_j.Time.else;
    Time.total           = Time.total           + Info_j.Time.total;  

    Log.param(j, :) = p_j(1);
    Log.cost(j, :) = [J_ocp_j, J_penalty_j];
    Log.KKT_error(j, :) = [KKT_error_primal_j, KKT_error_dual_j];
    if Option.recordLevel == 1
        Log.stepSize_primal(j, :) = [min(Info_j.Log.stepSize), sum(Info_j.Log.stepSize)/(Info_j.iterNum)];
        Log.stepSize_dual(j, :) = [min(Info_j.Log.stepSize), sum(Info_j.Log.stepSize)/(Info_j.iterNum)];
    end
    Log.VI_natural_residual(j) = VI_nat_res_j;
    Log.iterNum(j) = Info_j.iterNum;
    Log.timeElapsed(j) = Info_j.Time.total;
    if mod(j, 10) == 1
        disp('---------------------------------------------------------------------------------------------------------------------------------')
        headMsg = ' step  | param(mu)| cost(ocp/penalty) | KKT(primal/dual)| alpha_p(min/ave)| alpha_d(min/ave)| nat_res | iterNum | time(s) ';
        disp(headMsg)
    end
    prevIterMsg = [' ',...
        num2str(j,'%10.2d'), '/', num2str(continuationStepMaxNum,'%10.2d'),' |  ',...
        num2str(Log.param(j, 1), '%10.1e'), ' | ',...
        num2str(Log.cost(j, 1), '%10.2e'), ' ', num2str(Log.cost(j, 2), '%10.2e'),' | ',...
        num2str(Log.KKT_error(j, 1), '%10.1e'), ' ', num2str(Log.KKT_error(j, 2), '%10.1e'),' | ',...
        num2str(Log.stepSize_primal(j, 1), '%10.1e'), ' ', num2str(Log.stepSize_primal(j, 2), '%10.1e'), ' | ',...
        num2str(Log.stepSize_dual(j, 1), '%10.1e'), ' ' , num2str(Log.stepSize_dual(j, 2), '%10.1e'),' | ',...
        num2str(Log.VI_natural_residual(j), '%10.1e'),' |   ',...
        num2str(Log.iterNum(j), '%10.4d'),'  | ',...
        num2str(Log.timeElapsed(j), '%10.4f')];
    disp(prevIterMsg)

    %% step 3: check ternimation based on the current homotopy iterate
    solver_stats = Info_j.terminalStatus;
    if (solver_stats == 1) && (VI_nat_res_j <= VI_nat_res_tol)
        % IPOPT at this homotopy iteration finds the optimal solution satisfying the desired VI natural residual
        exitFlag = true;
        terminalStatus = 1;
        terminalMsg = Info_j.terminalMsg;
    elseif (solver_stats ~= 1)
        % IPOPT at this homotopy iteration fails to find the optimal solution
        exitFlag = true;
        terminalStatus = 0;
        terminalMsg = Info_j.terminalMsg;
    elseif j == continuationStepMaxNum
        % IPOPT still can not find the optimal solution in the final homotopy iteration
        exitFlag = true;
        terminalStatus = 0;
        terminalMsg = 'solver can not find the optimal solution satisfying the desired VI natural residual';
    else
        % IPOPT at this homotopy iteration (not the final) finds the optimal solution, prepare for next homotopy iteration
        exitFlag = false;
        % update initial guess
        z_Init_j = z_Opt_j;
        % update penalty parameter
        mu_j = p_j(1);
        mu_j = min([kappa_mu_times * mu_j, mu_End]);
        % update parameter vector
        p_j(1) = mu_j;
        % update continuation step counter
        j = j + 1;
    end
    %% step 4: check exitFlag and return optimal solution
    if exitFlag
        % return the current homotopy iterate as the optimal solution
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
        Info.time = sum(Log.timeElapsed);
        Info.iterNum = sum(Log.iterNum);
        Info.Log.param = Log.param(1 : j, :);
        Info.Log.cost = Log.cost(1 : j, :);
        Info.Log.KKT_error = Log.KKT_error(1 : j, :);
        Info.Log.stepSize_primal = Log.stepSize_primal(1 : j, :);
        Info.Log.stepSize_dual = Log.stepSize_dual(1 : j, :);
        Info.Log.VI_natural_residual = Log.VI_natural_residual(1 : j, :);
        Info.Log.iterNum = Log.iterNum(1 : j, :);
        Info.Log.timeElapsed = Log.timeElapsed(1 : j, :);
        % display homotopy terminal result and then break        
        disp('*--------------------------------------------- Solution Information ----------------------------------------------*')
        disp(['1. Terminal Message: ', Info.terminalMsg]) 
        disp('2. Continuation Step Message')
        disp(['- TimeElapsed: ................................. ', num2str(Info.time,'%10.3f'), ' s'])
        disp(['- Continuation Step: ........................... ', num2str(Info.continuationStepNum)])        
        disp(['- Time Per Continuation Step: .................. ', num2str(Info.time / Info.continuationStepNum,'%10.2f'), ' s/Step'])
        disp(['- Iterations: .................................. ', num2str(Info.iterNum)])
        disp(['- Time Per Iteration: .......................... ', num2str(1000 * Info.time / Info.iterNum,'%10.2f'), ' ms/Iter'])
        disp('3. Solution Message')
        disp(['- Cost(ocp): ................................... ', num2str(Info.cost.ocp,'%10.3e'), '; '])
        disp(['- Cost(penalty): ............................... ', num2str(Info.cost.penalty,'%10.3e'), '; '])
        disp(['- KKT(primal): ................................. ', num2str(Info.KKT_error.primal,'%10.3e'), '; '])
        disp(['- KKT(dual): ................................... ', num2str(Info.KKT_error.dual,'%10.3e')  '; '])
        disp(['- equilibrium constraint(natural residual): .... ', num2str(Info.VI_natural_residual,'%10.3e'), '; '])

        break
    end

end

end