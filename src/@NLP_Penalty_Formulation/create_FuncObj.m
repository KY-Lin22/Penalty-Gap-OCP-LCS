function FuncObj = create_FuncObj(self, nlp)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
FuncObj = struct();

%% NLP cost and constraint
% function
FuncObj.J = Function('J', {nlp.z, nlp.p}, {nlp.J}, {'z', 'p'}, {'J'});
FuncObj.h = Function('h', {nlp.z, nlp.p}, {nlp.h}, {'z', 'p'}, {'h'});
% jacobian
FuncObj.J_grad = Function('J_grad', {nlp.z, nlp.p}, {nlp.J_grad}, {'z', 'p'}, {'J_grad'});

%% other terms
% cost term
FuncObj.J_ocp = Function('J_ocp', {nlp.z, nlp.p}, {nlp.J_ocp}, {'z', 'p'}, {'J_ocp'});
FuncObj.J_penalty = Function('J_penalty', {nlp.z, nlp.p}, {nlp.J_penalty}, {'z', 'p'}, {'J_penalty'});
% D gap function, gradient and hessian
FuncObj.D_gap_func = Function('D_gap_func', {nlp.z}, {nlp.D_gap_func}, {'z'}, {'D_gap_func'});
FuncObj.D_gap_grad = Function('D_gap_grad', {nlp.z}, {nlp.D_gap_grad}, {'z'}, {'D_gap_grad'});
FuncObj.D_gap_hessian = Function('D_gap_hessian', {nlp.z}, {nlp.D_gap_hessian}, {'z'}, {'D_gap_hessian'});
end