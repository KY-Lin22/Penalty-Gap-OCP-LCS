function FuncObj = create_FuncObj(self, nlp)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
FuncObj = struct();

%% NLP cost and constraint
% cost and constraint function
FuncObj.J = Function('J', {nlp.z, nlp.p}, {nlp.J}, {'z', 'p'}, {'J'});
FuncObj.h = Function('h', {nlp.z, nlp.p}, {nlp.h}, {'z', 'p'}, {'h'});
FuncObj.c = Function('c', {nlp.z, nlp.p}, {nlp.c}, {'z', 'p'}, {'c'});

FuncObj.J_ocp = Function('J_ocp', {nlp.z, nlp.p}, {nlp.J_ocp}, {'z', 'p'}, {'J_ocp'});
FuncObj.J_penalty = Function('J_penalty', {nlp.z, nlp.p}, {nlp.J_penalty}, {'z', 'p'}, {'J_penalty'});

% cost and constraint jacobian
J_grad = jacobian(nlp.J, nlp.z);
FuncObj.J_grad = Function('J_grad', {nlp.z, nlp.p}, {J_grad}, {'z', 'p'}, {'J_grad'});

h_grad = jacobian(nlp.h, nlp.z);
h_grad_func = Function('h_grad_func', {nlp.z}, {h_grad}, {'z'}, {'h_grad'});
FuncObj.h_grad = sparse(h_grad_func(zeros(nlp.Dim.z, 1))); % constant matrix

c_grad = jacobian(nlp.c, nlp.z);
FuncObj.c_grad = Function('c_grad', {nlp.z, nlp.p}, {c_grad}, {'z', 'p'}, {'c_grad'});

J_ocp_grad = jacobian(nlp.J_ocp, nlp.z);
FuncObj.J_ocp_grad = Function('J_ocp_grad', {nlp.z, nlp.p}, {J_ocp_grad}, {'z', 'p'}, {'J_ocp_grad'});

J_penalty_grad = jacobian(nlp.J_penalty, nlp.z);
FuncObj.J_penalty_grad = Function('J_penalty_grad', {nlp.z, nlp.p}, {J_penalty_grad}, {'z', 'p'}, {'J_penalty_grad'});

% cost hessian
[J_ocp_hessian, ~] = hessian(nlp.J_ocp, nlp.z);
J_ocp_hessian_func = Function('J_ocp_hessian_func', {nlp.z}, {J_ocp_hessian}, {'z'}, {'J_ocp_hessian'});
FuncObj.J_ocp_hessian = sparse(J_ocp_hessian_func(zeros(nlp.Dim.z, 1))); % constant matrix

[J_penalty_hessian, ~] = hessian(nlp.J_penalty, nlp.z);
FuncObj.J_penalty_hessian = Function('J_penalty_hessian', {nlp.z, nlp.p}, {J_penalty_hessian}, {'z', 'p'}, {'J_penalty_hessian'});

% D gap function
if isfield(nlp, 'D_gap_func')
    FuncObj.D_gap_func = Function('D_gap_func', {nlp.z}, {nlp.D_gap_func}, {'z'}, {'D_gap_func'});
end

end