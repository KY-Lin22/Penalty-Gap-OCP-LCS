clear all
clc

%% create OCP, NLP, and solver
% construct OCP problem
OCP = OCP_Vieira_LCS_analytic();
Option.penalty_problem = 'complementarity_based'; % 'gap_based', 'complementarity_based'
Option.gap_penalty_strategy = 'func_direct'; % 'func_direct', 'grad_cvx_over_nonlinear'
Option.CHKS_param = 1e-5;
Option.Huber_param = 0.1;
Option.D_gap_param_a = 0.5;
Option.D_gap_param_b = 1.5;
NLP = NLP_Penalty_Formulation(OCP, Option);

solver = Penalty_IPOPT_Based_Solver(OCP, NLP);

%% set option and solve problem
% set option
solver.Option.Homotopy.kappa_mu_times = 1.2;
solver.Option.Homotopy.VI_nat_res_tol = 1e-2;
%
z_Init = ones(NLP.Dim.z, 1);
p_Init = 1;
p_End = 5;
[z_Opt, Info] = solver.solve_NLP(z_Init, p_Init, p_End);

%% show result
plotResult_Vieira_LCS_analytic(OCP, NLP, z_Opt)