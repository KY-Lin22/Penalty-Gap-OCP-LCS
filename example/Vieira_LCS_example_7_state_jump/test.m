clear all
clc

%% create OCP, NLP, and solver
% construct OCP problem
OCP = OCP_Vieira_LCS_state_jump();
% construct NLP problem
Option.reformulation_strategy = 'penalty'; % 'relaxation', 'penalty', 'smoothing'
Option.relaxation_problem = 'Lin_Fukushima'; % 'Scholtes', 'Lin_Fukushima'
Option.penalty_problem = 'gap_based'; % 'gap_based', 'complementarity_based'
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

z_Init = ones(NLP.Dim.z, 1);
% p = 1e1;
% [z_Opt, Info] = solver.solve_NLP_single(z_Init, p);

p_Init = 1e1;
p_End = 1e5;
num_test = 5;
time_rec = 0;
for i = 1 : num_test
    [z_Opt, Info] = solver.solve_NLP(z_Init, p_Init, p_End);
    time_rec = time_rec + Info.time;
end
disp(['average time: ', num2str(time_rec/num_test)])
disp(['cost: ', num2str(Info.cost.ocp)])

%% show result
plotResult_Vieira_LCS_state_jump(OCP, NLP, z_Opt)
