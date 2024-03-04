function polar_func = create_polar_func(self, OCP)
%UNTITLED24 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
% load parameter
param_b = self.D_gap_param_b;
% symbolic variable
theta = SX.sym('theta', OCP.Dim.lambda, 1);
r = SX.sym('r', OCP.Dim.lambda, 1);
% polar transfer
lambda = r .* cos(theta + atan(param_b));
eta = r .* sin(theta + atan(param_b));

polar_func = Function('polar_func', {theta, r}, {lambda, eta}, {'theta', 'r'}, {'lambda', 'eta'});
end