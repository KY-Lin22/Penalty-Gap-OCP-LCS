clear all
clc

%% create OCP, NLP, and solver
% construct OCP problem
OCP = OCP_Vieira_LCS_lambda_penalty();
% construct NLP problemLin_Fukushima
Option.reformulation_strategy = 'smoothing'; % 'relaxation', 'penalty', 'smoothing'
Option.relaxation_problem = 'Lin_Fukushima'; % 'Scholtes', 'Lin_Fukushima'
Option.penalty_problem = 'gap_based'; % 'gap_based', 'complementarity_based'
Option.CHKS_param = 0;
Option.D_gap_param_a = 0.9;
Option.D_gap_param_b = 1.1;
NLP = NLP_Formulation(OCP, Option);
% construct solver
switch NLP.reformulation_strategy
    case 'relaxation'
        solver = IPOPT_Based_Solver(OCP, NLP);
    case 'penalty'
        switch NLP.penalty_problem
            case 'gap_based'
                solver = Penalty_Gap_Solver(OCP, NLP);
            case 'complementarity_based'
                solver = IPOPT_Based_Solver(OCP, NLP);
        end
    case 'smoothing'
        solver = IPOPT_Based_Solver(OCP, NLP);
end

%% set option and solve problem
solver.Option.Homotopy.kappa_mu_times = 1.2;
solver.Option.Homotopy.VI_nat_res_tol = 1e-2;

z_Init = randn(NLP.Dim.z, 1);
% p = 1e1;
% [z_Opt, Info] = solver.solve_NLP_single(z_Init, p);

p_Init = 1e1;
p_End = 1e5;
[z_Opt, Info] = solver.solve_NLP(z_Init, p_Init, p_End);

%% show result
plotResult_Vieira_LCS_lambda_penalty(OCP, NLP, z_Opt)
