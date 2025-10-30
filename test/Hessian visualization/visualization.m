clear all
clc
import casadi.*

% parameter
param_a = 0.5;
param_b = 2;
dim_lambda = 5;
% variable
lambda = SX.sym('lambda', dim_lambda, 1);
eta = SX.sym('eta', dim_lambda, 1);
xi_vec = [lambda; eta];
xi_ele = reshape([lambda'; eta'], 2 * dim_lambda, 1);
% D-gap func
phi_ab = (1/2*param_a) * (norm(eta, 2)^2 - norm(max(zeros(dim_lambda, 1), eta - param_a * lambda), 2)^2) ...
    - (1/2*param_b) * (norm(eta, 2)^2 - norm(max(zeros(dim_lambda, 1), eta - param_b * lambda), 2)^2);
% hessian
[phi_ab_hess_vec, ~] = hessian(phi_ab, xi_vec);
[phi_ab_hess_ele, ~] = hessian(phi_ab, xi_ele);

figure(1)
spy(phi_ab_hess_vec, 30)
%title('Hessian (vector)')
figure(2)
spy(phi_ab_hess_ele, 30)
%title('Hessian (element)')