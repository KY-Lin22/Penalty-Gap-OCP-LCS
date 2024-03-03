function regular_func = create_regular_func(self)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
% load parameter
param_a = self.D_gap_param_a;
param_b = self.D_gap_param_b;
% symbolic variable
u = SX.sym('u', 1, 1);
v = SX.sym('v', 1, 1);
f = param_a/2*u^2 - u*v + 1/(2*param_b)*v^2; % negative D gap concave part
regular_func = Function('regular_func', {u, v}, {f}, {'u', 'v'}, {'f'});
end