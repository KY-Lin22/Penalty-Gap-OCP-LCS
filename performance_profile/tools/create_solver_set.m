function solver_set = create_solver_set(OCP_func_handle, nStages_sequ, NLP_option_set, solver_Option)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%% generate OCP problem set
OCP_problem_set = cell(numel(OCP_func_handle), numel(nStages_sequ));
for i = 1 : numel(OCP_func_handle)
    for j = 1 : numel(nStages_sequ)
        OCP_i_j = OCP_func_handle{i}();
        OCP_i_j.nStages = nStages_sequ{j};
        OCP_i_j.timeStep = OCP_i_j.TimeHorizon ./ OCP_i_j.nStages;
        OCP_problem_set{i, j} = OCP_i_j;
    end
end
OCP_problem_set = reshape(OCP_problem_set, [], 1);

%% generate NLP problem set
NLP_problem_set = cell(numel(OCP_problem_set), numel(NLP_option_set));
for i = 1 : numel(OCP_problem_set)
    for j = 1 : numel(NLP_option_set)
        OCP_i = OCP_problem_set{i};
        NLP_option_j = NLP_option_set{j};
        NLP_i_j = NLP_Formulation(OCP_i, NLP_option_j);
        NLP_problem_set{i, j} = NLP_i_j;
    end
end

%% generate solver set
solver_set = cell(size(NLP_problem_set));
for i = 1 : size(NLP_problem_set, 1)
    for j = 1 : size(NLP_problem_set, 2)
        OCP_i = OCP_problem_set{i};
        NLP_i_j = NLP_problem_set{i, j};
        switch NLP_i_j.reformulation_strategy
            case 'relaxation'
                solver_i_j = IPOPT_Based_Solver(OCP_i, NLP_i_j);
            case 'penalty'
                switch NLP_i_j.penalty_problem
                    case 'gap_based'
                        solver_i_j = Penalty_Gap_Solver(OCP_i, NLP_i_j);
                    case 'complementarity_based'
                        solver_i_j = IPOPT_Based_Solver(OCP_i, NLP_i_j);
                end
            case 'smoothing'
                solver_i_j = IPOPT_Based_Solver(OCP_i, NLP_i_j);
        end
        % homotopy option
        solver_i_j.Option.Homotopy.kappa_mu_times = solver_Option.kappa_mu_times;
        solver_i_j.Option.Homotopy.VI_nat_res_tol = solver_Option.VI_nat_res_tol;
        solver_set{i, j} = solver_i_j;
    end
end

end