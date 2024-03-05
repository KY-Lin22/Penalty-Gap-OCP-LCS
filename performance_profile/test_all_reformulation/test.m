close all
clear all
clc

% run test (it needs about )
% [Rec, NLP_reformulation_name] = run_test_all_reformulation();

save('Data_performance_test.mat', 'Rec', 'NLP_reformulation_name')

%%
Data_performance_test = load('Data_performance_test.mat');

%
performance_matrix = Data_performance_test.Rec.time;
solver_name = Data_performance_test.NLP_reformulation_name;

plot_performance_profile(performance_matrix(:, 1:4), solver_name(1:4))

