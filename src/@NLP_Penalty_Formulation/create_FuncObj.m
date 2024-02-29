function FuncObj = create_FuncObj(self, nlp)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
FuncObj = struct();

%% function
% cost J
FuncObj.J = Function('J', {nlp.z, nlp.p}, {nlp.J}, {'z', 'p'}, {'J'});
% constraint h
FuncObj.h = Function('h', {nlp.z, nlp.p}, {nlp.h}, {'z', 'p'}, {'h'});
% cost term
FuncObj.J_ocp = Function('J_ocp', {nlp.z, nlp.p}, {nlp.J_ocp}, {'z', 'p'}, {'J_ocp'});
FuncObj.J_penalty = Function('J_penalty', {nlp.z, nlp.p}, {nlp.J_penalty}, {'z', 'p'}, {'J_penalty'});
% D gap function
FuncObj.D_gap = Function('D_gap', {nlp.z}, {nlp.D_gap}, {'z'}, {'D_gap'});

%% Jacobian and Hessian
% cost gradient
FuncObj.J_grad = Function('J_grad', {nlp.z, nlp.p}, {nlp.J_grad}, {'z', 'p'}, {'J_grad'});
% D_gap_grad
FuncObj.D_gap_grad = Function('D_gap_grad', {nlp.z}, {nlp.D_gap_grad}, {'z'}, {'D_gap_grad'});

end