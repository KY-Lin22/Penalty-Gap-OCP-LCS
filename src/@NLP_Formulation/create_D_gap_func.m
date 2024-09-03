function D_gap_func = create_D_gap_func(self)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
% load parameter
param_a = self.D_gap_param_a;
param_b = self.D_gap_param_b;
% symbolic variable
u = SX.sym('u', 1, 1);
v = SX.sym('v', 1, 1);
% D gap function
phi_ab = (param_b - param_a)/(2*param_a*param_b)*v^2 ...
    - 1/(2*param_a)*(max(0, v - param_a*u))^2 ...
    + 1/(2*param_b)*(max(0, v - param_b*u))^2;
% output
D_gap_func = Function('D_gap_func', {u, v}, {phi_ab}, {'u', 'v'}, {'phi_ab'});

end