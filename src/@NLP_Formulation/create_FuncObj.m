function FuncObj = create_FuncObj(self, nlp)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
FuncObj = struct();

%% NLP function
% cost and constraint function
FuncObj.J = Function('J', {nlp.z, nlp.p}, {nlp.J}, {'z', 'p'}, {'J'});
FuncObj.h = Function('h', {nlp.z, nlp.p}, {nlp.h}, {'z', 'p'}, {'h'});
FuncObj.c = Function('c', {nlp.z, nlp.p}, {nlp.c}, {'z', 'p'}, {'c'});

FuncObj.J_ocp = Function('J_ocp', {nlp.z, nlp.p}, {nlp.J_ocp}, {'z', 'p'}, {'J_ocp'});
FuncObj.J_penalty = Function('J_penalty', {nlp.z, nlp.p}, {nlp.J_penalty}, {'z', 'p'}, {'J_penalty'});

%% Function object for penalty gap solver
if strcmp(self.reformulation_strategy, 'penalty') && strcmp(self.penalty_problem, 'gap_based')
    % flag vector
    FuncObj.w = Function('w', {nlp.z}, {nlp.w_formula}, {'z'}, {'w'});
    % NLP Jacobian
    J_grad = jacobian(nlp.J, nlp.z);
    h_grad = jacobian(nlp.h, nlp.z);
    % NLP Hessian    
    [J_hessian, ~] = hessian(nlp.J, nlp.z);
    [regular_hessian, ~] = hessian(nlp.regular_func_w, nlp.z);
    J_hessian_regular = J_hessian + regular_hessian;
    % Lagrangian function and jacobian
    gamma_h = SX.sym('gamma_h', nlp.Dim.h, 1);
    LAG = nlp.J + gamma_h' * nlp.h;
    LAG_grad = jacobian(LAG, nlp.z);
    FuncObj.LAG_grad = Function('LAG_grad', {nlp.z, gamma_h, nlp.p}, {LAG_grad}, {'z', 'gamma_h', 'p'}, {'LAG_grad'});
    % KKT residual
    KKT_residual = [J_grad'; nlp.h];
    FuncObj.KKT_residual = Function('KKT_residual', {nlp.z, nlp.p}, {KKT_residual}, {'z', 'p'}, {'KKT_residual'});
    % KKT matrix
    KKT_matrix = [J_hessian_regular, h_grad';...
                h_grad, SX(nlp.Dim.h, nlp.Dim.h)];
    % Newton step
    dz_gamma_h = -KKT_matrix\KKT_residual;
    FuncObj.dz_gamma_h = Function('dz_gamma_h', {nlp.z, nlp.p, nlp.w}, {dz_gamma_h}, {'z', 'p', 'w'}, {'dz_gamma_h'});
end

end