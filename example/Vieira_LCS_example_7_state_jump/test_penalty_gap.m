clear all
clc

%% create OCP, NLP, and solver
% construct OCP problem
OCP = OCP_Vieira_LCS_state_jump();
Option.penalty_problem = 'gap_based'; % 'gap_based', 'complementarity_based'
Option.CHKS_param = 0;
Option.D_gap_param_a = 0.1;
Option.D_gap_param_b = 10;
NLP = NLP_Penalty_Formulation(OCP, Option);

solver = Penalty_Gap_Solver(OCP, NLP);

%% set option and solve problem
solver.Option.maxIterNum = 500;
solver.Option.tol.KKT_error_primal = 1e-6; 
solver.Option.tol.KKT_error_dual = 1e-4; 
solver.Option.tol.KKT_error_total = 1e-6; 
solver.Option.tol.dzNorm = 1e-8;
solver.Option.penalty_hessian_regularization = 1;
solver.Option.Homotopy.kappa_mu_times = 1.5;
solver.Option.Homotopy.VI_nat_res_tol = 1e-2;
z_Init = randn(NLP.Dim.z, 1);
% p = 0.1;
% [z_Opt, Info] = solver.solve_NLP_single(z_Init, p);
p_Init = 0.01;
p_End = 100;
[z_Opt, Info] = solver.solve_NLP(z_Init, p_Init, p_End);

%% show result
plotResult_Vieira_LCS_state_jump(OCP, NLP, z_Opt)
