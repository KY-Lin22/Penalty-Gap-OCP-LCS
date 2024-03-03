clear all
clc

%% create OCP, NLP, and solver
% construct OCP problem
OCP = OCP_Vieira_LCS_analytic();
Option.penalty_problem = 'complementarity_based'; % 'gap_based', 'complementarity_based'
Option.CHKS_param = 0;
Option.D_gap_param_a = 0.95;
Option.D_gap_param_b = 1;
NLP = NLP_Penalty_Formulation(OCP, Option);

solver = Penalty_IPOPT_Based_Solver(OCP, NLP);

%% set option and solve problem
% set option
solver.Option.Homotopy.kappa_mu_times = 1.2;
solver.Option.Homotopy.VI_nat_res_tol = 1e-2;
%
z_Init = ones(NLP.Dim.z, 1);
p_Init = 10;
p_End = 20;
[z_Opt, Info] = solver.solve_NLP(z_Init, p_Init, p_End);

%% show result
plotResult_Vieira_LCS_analytic(OCP, NLP, z_Opt)