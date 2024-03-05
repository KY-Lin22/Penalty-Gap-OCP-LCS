function Rec = run_solver_test(solver_set, p_Init, p_End)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% record cost and time
Rec.cost = zeros(size(solver_set));
Rec.time = zeros(size(solver_set));
testNum = 5;
ignoreNum = 2;
% solve
for i = 1 : size(solver_set, 1)
    for j = 1 : size(solver_set, 2)
        disp(['-----------------------------',...
            'prob: ', num2str(i), ' / ', num2str(size(solver_set, 1)), '; ', ...
            'solver: ',  num2str(j), ' / ', num2str(size(solver_set, 2)),...
            '-----------------------------'])
        solver_i_j = solver_set{i, j};
        z_Init_i_j = ones(solver_i_j.NLP.Dim.z, 1);
        time_rec = zeros(1, testNum + ignoreNum);
        for k = 1 : testNum + ignoreNum
            [~, Info_i_j] = solver_i_j.solve_NLP(z_Init_i_j, p_Init, p_End);
            time_rec(1, k) = Info_i_j.time;
        end      
        if Info_i_j.terminalStatus == 1
            Rec.cost(i, j) = Info_i_j.cost.ocp;
            Rec.time(i, j) = sum(time_rec(1, ignoreNum + 1 : end))/testNum;
        else
            % set fail case as inf
            Rec.cost(i, j) = inf;
            Rec.time(i, j) = inf;
        end
    end
end

end