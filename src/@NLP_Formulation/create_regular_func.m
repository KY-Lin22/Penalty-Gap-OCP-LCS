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
% 1. using negative D gap concave part (give up)
% f = param_a/2*u^2 - u*v + 1/(2*param_b)*v^2; 
% 2. using negative eigenvalue to construct a quadratic term
nega_eig = min(eig([-param_a, 1; 1, -1/param_b]));
regular_param = -nega_eig;
f = 0.5*regular_param*(u^2 + v^2);
regular_func = Function('regular_func', {u, v}, {f}, {'u', 'v'}, {'f'});
end