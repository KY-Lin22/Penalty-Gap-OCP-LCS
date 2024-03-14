function [group_func, decouple_func_lambda, decouple_func_eta] = create_element_wise_concatenation_func(self, OCP)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
import casadi.*

%% group function
u_group = SX.sym('u_group', OCP.Dim.lambda, 1);
v_group = SX.sym('v_group', OCP.Dim.lambda, 1);

xi_group = reshape([u_group'; v_group'], 2 * OCP.Dim.lambda, 1);
group_func = Function('group_func', {u_group, v_group}, {xi_group}, ...
    {'u_group', 'v_group'}, {'xi_group'});

%% decouple function
xi_decouple = SX.sym('xi_decouple', 2 * OCP.Dim.lambda, 1);
xi_decouple_reshaped = reshape(xi_decouple, 2, OCP.Dim.lambda);
u_decouple = xi_decouple_reshaped(1, :)';
v_decouple = xi_decouple_reshaped(2, :)';

decouple_func_lambda = Function('decouple_func_lambda', {xi_decouple}, {u_decouple}, ...
    {'xi_decouple'}, {'u_decouple'});
decouple_func_eta = Function('decouple_func_eta', {xi_decouple}, {v_decouple}, ...
    {'xi_decouple'}, {'v_decouple'});
end