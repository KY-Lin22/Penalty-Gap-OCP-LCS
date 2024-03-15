close all
clear all
clc

% run test (it needs about )
% [Rec, NLP_reformulation_name] = run_test_all_reformulation();

save('Data_performance_test.mat', 'Rec', 'NLP_reformulation_name')

%%
Data_performance_test = load('Data_performance_test.mat');

%
performance_matrix = Data_performance_test.Rec.cost;
solver_name = Data_performance_test.NLP_reformulation_name;

plot_performance_profile(performance_matrix(:, 2:6), solver_name(2:6))

