function [D_gap_func, D_gap_grad] = create_D_gap_func(self)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
% load parameter
CHKS_param = self.CHKS_param;
param_a = self.D_gap_param_a;
param_b = self.D_gap_param_b;
% symbolic variable
u = SX.sym('u', 1, 1);
v = SX.sym('v', 1, 1);
% smoothing max(0, v-bu) and max(0, v-au) using 0.5*(sqrt(x^2 + 4*CHKS_param^2) + x) for max(0, x) 
smooth_max_0_vau = 0.5*(sqrt((v - param_a*u)^2 + 4*CHKS_param^2) + (v - param_a*u));
smooth_max_0_vbu = 0.5*(sqrt((v - param_b*u)^2 + 4*CHKS_param^2) + (v - param_b*u));
% D gap function
phi_ab = (param_b - param_a)/(2*param_a*param_b)*v^2 ...
    - 1/(2*param_a)*(smooth_max_0_vau)^2 ...
    + 1/(2*param_b)*(smooth_max_0_vbu)^2;
% D gap function gradient
phi_ab_grad_u = smooth_max_0_vau - smooth_max_0_vbu;
phi_ab_grad_v = (param_b - param_a)/(param_a*param_b)*v ...
    - 1/(param_a)*(smooth_max_0_vau) ...
    + 1/(param_b)*(smooth_max_0_vbu);
phi_ab_grad = [phi_ab_grad_u, phi_ab_grad_v];

% output
D_gap_func = Function('D_gap_func', {u, v}, {phi_ab}, {'u', 'v'}, {'phi_ab'});
D_gap_grad = Function('D_gap_grad', {u, v}, {phi_ab_grad}, {'u', 'v'}, {'phi_ab_grad'});

end