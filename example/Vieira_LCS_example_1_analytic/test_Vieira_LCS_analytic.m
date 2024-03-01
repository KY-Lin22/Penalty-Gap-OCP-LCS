clear all
clc

%% create OCP, NLP, and solver
% construct OCP problem
OCP = OCP_Vieira_LCS_analytic();
Option.CHKS_param = 1e-5;
Option.Huber_param = 0.1;
Option.D_gap_param_a = 0.9;
Option.D_gap_param_b = 1.1;
NLP = NLP_Penalty_Formulation(OCP, Option);

solver = Penalty_Gap_Solver(OCP, NLP);

%% set option and solve problem
% one solve
solver.Option.maxIterNum = 500;
solver.Option.tol.KKT_error_primal = 1e-6; 
solver.Option.tol.KKT_error_dual = 1e-4; 
solver.Option.tol.KKT_error_total = 1e-6; 
solver.Option.tol.dzNorm = 1e-8;
solver.Option.recordLevel = 1;
solver.Option.printLevel = 2;
solver.Option.qpSolver = 'osqp'; % 'osqp'
solver.Option.Homotopy.kappa_mu_times = 1.2;
solver.Option.Homotopy.VI_nat_res_tol = 1e-2;
z_Init = ones(NLP.Dim.z, 1);
% p = 1;
% [z_Opt, Info] = solver.solve_NLP_single(z_Init, p);
p_Init = 1;
p_End = 5;
[z_Opt, Info] = solver.solve_NLP(z_Init, p_Init, p_End);

%% show result
plotResult_Vieira_LCS_analytic(OCP, NLP, z_Opt)
