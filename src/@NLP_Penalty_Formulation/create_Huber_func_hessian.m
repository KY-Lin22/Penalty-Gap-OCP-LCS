function [Huber_func, Huber_hessian] = create_Huber_func_hessian(self, OCP)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*

% create pseudo Huber loss function
v = SX.sym('v', 1, 2 * OCP.Dim.lambda * OCP.nStages);
Huber_param = self.Huber_param;
Huber_formula = sum(sqrt(v.^2 + Huber_param^2) - Huber_param);
[Huber_hessian_formula, ~] = hessian(Huber_formula, v);

Huber_func = Function('Huber_func', {v}, {Huber_formula}, {'v'}, {'Huber_formula'});
Huber_hessian = Function('Huber_hessian', {v}, {Huber_hessian_formula}, {'v'}, {'Huber_hessian_formula'});
end