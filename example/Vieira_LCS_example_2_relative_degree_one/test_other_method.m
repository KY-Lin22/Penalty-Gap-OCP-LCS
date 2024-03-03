clear all
clc

%% create OCP, NLP, and solver
% construct OCP problem
OCP = OCP_Vieira_LCS_relative_degree_one();
Option.penalty_problem = 'gap_based'; % 'gap_based', 'complementarity_based'
Option.CHKS_param = 0;
Option.D_gap_param_a = 0.95;
Option.D_gap_param_b = 1;
NLP = NLP_Penalty_Formulation(OCP, Option);

solver = Penalty_IPOPT_Based_Solver(OCP, NLP);

%% set option and solve problem
% set option
solver.Option.Homotopy.kappa_mu_times = 1.5;
solver.Option.Homotopy.VI_nat_res_tol = 1e-2;
%
z_Init = randn(NLP.Dim.z, 1);
p_Init = 1;
p_End = 10;
[z_Opt, Info] = solver.solve_NLP(z_Init, p_Init, p_End);

%% show result
plotResult_Vieira_LCS_relative_degree_one(OCP, NLP, z_Opt)