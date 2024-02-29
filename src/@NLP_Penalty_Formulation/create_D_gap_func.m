function D_gap_func = create_D_gap_func(self, OCP)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
import casadi.*

%% convex function d and its jacobian
x = SX.sym('x', OCP.Dim.lambda, 1); 
d = 0.5 * x' * eye(OCP.Dim.lambda) * x;
dx = jacobian(d, x);
d_func = Function('d_func', {x}, {d}, {'x'}, {'d'});
d_grad = Function('d_grad', {x}, {dx}, {'x'}, {'dx'});

%% weighting generalized primal gap function phi_c
lambda = SX.sym('lambda', OCP.Dim.lambda, 1);
eta = SX.sym('eta', OCP.Dim.lambda, 1); % auxiliary variable for complementarity function in gap function

c = SX.sym('c', 1, 1); % weighting nonnegative scalar parameter 
CHKS_param = self.CHKS_param;

stationary_point_c = lambda - (1/c) * eta;
omega_c = 0.5*(sqrt(stationary_point_c.^2 + 4 * CHKS_param^2) + stationary_point_c);
p_c = d_func(lambda) - d_func(omega_c) + d_grad(lambda) * (omega_c - lambda);
phi_c = eta' * (lambda - omega_c) + c * p_c;
phi_c_func = Function('phi_c_func', {lambda, eta, c}, {phi_c}, {'lambda', 'eta', 'c'}, {'phi_c'});

%% D gap function
param_a = self.D_gap_param_a;
param_b = self.D_gap_param_b;
phi_a = phi_c_func(lambda, eta, param_a);
phi_b = phi_c_func(lambda, eta, param_b);
phi_ab = phi_a - phi_b;

% output
D_gap_func = Function('D_gap_func', {lambda, eta}, {phi_ab}, {'lambda', 'eta'}, {'phi_ab'});

end