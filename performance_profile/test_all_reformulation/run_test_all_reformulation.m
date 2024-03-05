function [Rec, NLP_reformulation_name] = run_test_all_reformulation()
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%% OCP example to be tested 
OCP_func_handle = {...
    @Vieira_LCS_Analytic_1;...
    @Vieira_LCS_Analytic_2;...
    @Vieira_LCS_Rel_Deg_One;...
    @Vieira_LCS_High_Dim;...
    @Vieira_LCS_Without_Penalty;...
    @Vieira_LCS_With_Penalty_1;...
    @Vieira_LCS_With_Penalty_2;...
    @Vieira_LCS_With_Penalty_3;...
    @Vieira_LCS_Control_Jump;...
    @Vieira_LCS_State_Jump_1;...
    @Vieira_LCS_State_Jump_2};

nStages_sequ = {50, 80, 100, 200, 250, 400};

%% penalty NLP reformulation to be tested
NLP_option_set = {...
    struct('penalty_problem', 'complementarity_based'),...
    struct('penalty_problem', 'gap_based',...
    'CHKS_param', 0, 'D_gap_param_a', 0.9, 'D_gap_param_b', 1.1),...
    struct('penalty_problem', 'gap_based',...
    'CHKS_param', 1e-5, 'D_gap_param_a', 0.9, 'D_gap_param_b', 1.1),...
    };
NLP_reformulation_name = {'complementarity', 'Gap (without smooth)', 'Gap (with smooth)'};

%% solver option
solver_Option.kappa_mu_times = 1.2;
solver_Option.VI_nat_res_tol = 1e-2;

%% generate solver set
solver_set = create_solver_set(OCP_func_handle, nStages_sequ, NLP_option_set, solver_Option);

%% run test
% parameter
mu_Init = 1e1;
mu_End = 1e5;
p_Init = mu_Init;
p_End = mu_End;
% solve
Rec = run_solver_test(solver_set, p_Init, p_End);

end