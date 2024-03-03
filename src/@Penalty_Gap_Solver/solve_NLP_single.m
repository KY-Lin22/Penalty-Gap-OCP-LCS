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
OCP = self.OCP;
NLP = self.NLP;
Option = self.Option;
% check input z_Init and p
if ~all(size(z_Init) == [NLP.Dim.z, 1])
    error('z_Init has wrong dimension')
end
if ~all(size(p) == [NLP.Dim.p, 1])
    error('p has wrong dimension')
end
% check parameter
if (p(1) < 0)
    error('penalty parameter mu (i.e, p_1) should be nonnegative')
end

%% iteration routine (z: previous iterate z_{k-1}, z_k: current iterate z_{k}) 
% time record
Time = struct('gradEval', 0, 'searchDirection', 0, 'else', 0, 'total', 0);

k = 1;
z = z_Init;
gamma_h = zeros(NLP.Dim.h, 1);
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

    % constraint
    h = full(NLP.FuncObj.h(z, p));
    % cost Jacobian
    J_grad = full(NLP.FuncObj.J_grad(z, p));
    % penalty hessian (exact)
    J_penalty_hessian = sparse(NLP.FuncObj.J_penalty_hessian(z, p));
    % penalty hessian (regularization)
    Z = reshape(z, NLP.Dim.z_Node(end), OCP.nStages);
    LAMBDA = Z(NLP.Dim.z_Node(2) + 1 : NLP.Dim.z_Node(3), :);
    ETA = Z(NLP.Dim.z_Node(3) + 1 : NLP.Dim.z_Node(4), :);
    lambda = reshape(LAMBDA, [], 1);
    eta = reshape(ETA, [], 1); 
    w = ((NLP.D_gap_param_b*lambda > eta) .* (eta > NLP.D_gap_param_a*lambda))';
    regular_hessian = sparse(NLP.FuncObj.regular_hessian(z, p, w));
    J_penalty_hessian = J_penalty_hessian + regular_hessian;
    % KKT error (L_inf norm)
    KKT_error_primal = norm(h, inf);
    KKT_error_dual = norm(J_grad + gamma_h' * NLP.FuncObj.h_grad, inf);

    timeElasped_gradEval = toc(timeStart_gradEval);

    %% step 2: search direction evaluation based on previous iterate z
    timeStart_searchDirection = tic;

    % KKT residual
    KKT_residual = [-J_grad'; -h];
    % KKT matrix
    [i_J_pen_hess, j_J_pen_hess, s_J_pen_hess] = find(J_penalty_hessian);
    KKT_matrix_update = sparse(i_J_pen_hess, j_J_pen_hess, s_J_pen_hess,...
        NLP.Dim.z + NLP.Dim.h, NLP.Dim.z + NLP.Dim.h, length(s_J_pen_hess));
    KKT_matrix = self.KKT_matrix_constant + KKT_matrix_update;
    % solve linear system
    dz_gamma_h_k = KKT_matrix\KKT_residual;
    % extract primal and dual part
    dz = dz_gamma_h_k(1 : NLP.Dim.z, 1);
    gamma_h_k = dz_gamma_h_k(NLP.Dim.z + 1 : end, 1);
    dzNorm = norm(dz, inf);

    timeElasped_searchDirection = toc(timeStart_searchDirection);

    %% step 3: check whether we can terminate successfully based on the previous iterate z
    if max([KKT_error_primal, KKT_error_dual]) < Option.tol.KKT_error_total
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

    %% step 4: record and print information of this iteration k
    timeElasped_total = toc(timeStart_total);

    Time.gradEval = Time.gradEval + timeElasped_gradEval;
    Time.searchDirection = Time.searchDirection + timeElasped_searchDirection;   
    Time.total = Time.total + timeElasped_total;
    % print
    if Option.printLevel == 2
        % head
        if mod(k, 10) == 1
            disp('----------------------------------------------------------------------------------------')
            headMsg = ' Iter | cost(ocp/penalty) |  gapRes  | KKT(primal/dual)|  dzNorm  | time(ms) |';
            disp(headMsg)
        end
        % previous iterate message
        prevIterMsg = [' ',...
            num2str(k,'%10.3d'),'  | ',...
            num2str(full(NLP.FuncObj.J_ocp(z, p)), '%10.2e'), ' ',...
            num2str(full(NLP.FuncObj.J_penalty(z, p)), '%10.2e'),' | ',...
            num2str(norm(full(NLP.FuncObj.D_gap_func(z)), inf), '%10.2e'), ' | ',...
            num2str(KKT_error_primal, '%10.1e'), ' ',...
            num2str(KKT_error_dual, '%10.1e'),' | ',...
            num2str(dzNorm,'%10.2e'), ' | ',...
            num2str(1000 * timeElasped_total,'%10.2e'), ' | '];
        disp(prevIterMsg)
    end

    %% step 5: prepare next iteration
    k = k + 1;
    z = z + dz;
    gamma_h = gamma_h_k;
end

%% return optimal solution and create information
% return previous iterate as solution
z_Opt = z;
% create Info (basic: time, iterNum, terminalStatus)
Time.else = Time.total - Time.searchDirection - Time.gradEval;
Info.Time = Time;
Info.iterNum = k - 1;
Info.terminalStatus = terminalStatus;
Info.terminalMsg = terminalMsg;
% create Info (corresponds to the solution z_Opt: dual variable, cost, KKT, natural residual)
Info.gamma_h             = gamma_h;
Info.cost_ocp            = full(NLP.FuncObj.J_ocp(z, p));
Info.cost_penalty        = full(NLP.FuncObj.J_penalty(z, p));
Info.KKT_error_primal    = KKT_error_primal;
Info.KKT_error_dual      = KKT_error_dual;
Info.VI_natural_residual = self.evaluate_natural_residual(z);
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
    disp(['- KKT (primal): ..............', num2str(Info.KKT_error_primal,'%10.3e'), '; '])
    disp(['- KKT (dual): ................', num2str(Info.KKT_error_dual,'%10.3e'), '; '])
    disp(['- VI natural residual: .......', num2str(Info.VI_natural_residual,'%10.3e'), '; '])
end

end