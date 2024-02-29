function nlp = create_penalty_gap_based_NLP(self, OCP)
% formulate a penalty gap function based NLP
%
% OCP_LCS has the form:
%  min  int_0^T 0.5 * (x^T Q_x x + u^T Q_u u + lambda^T Q_lambda lambda) dt,
%  s.t. Dot{x} = A x + B u + E lambda
%       eta = C x + D u + F lambda
%       0 <= lambda \perp eta >= 0
%
% NLP has the form:
%  min  J(z, p),
%  s.t. h(z, p) = 0,
% where: (1) z: collects all the variables to be optimized and arranged in a stagewise manner
%            z = [z_1;...z_n;...z_N] and z_n = [x_n; u_n; lambda_n; eta_n] 
%            with x_n:      system state
%                 u_n:      control                 
%                 lambda_n: algebraic variable   
%                 eta_n:  auxiliary variable for complementarity function
%        (2) p: collects all the NLP problem parameters p = [mu]
%            with mu: nonnegative penalty parameter for penalty function
%        (3) J: cost function J = J_ocp + J_penalty
%            J_ocp = sum(J_stage_n)*dt
%            J_penalty = mu*huber_func([D_gap_1, ... D_gap_n, ...D_gap_N])
%            with J_stage_n:   stage cost defined in OCP
%        (4) h: equality constraint arranged in a stagewise manner
%            h = [h_1;...h_n;...h_N] and 
%            h_n = [x_{n-1} - x_n + f(x_n, u_n, lambda_n)*dt;
%                   g(x_n, u_n, lambda_n) - eta_n]     
%
% output: nlp is a structure with fields:
%         z: variable
%         p: parameter
%         J, h: cost and constraint function,  
%         J_ocp, J_penalty: cost term
%         J_grad: cost term gradient
%         h_grad: constraint Jacobian (constant)
%         J_ocp_hessian: ocp cost term hessian (constant)
%         D_gap_grad: used to evaluate J_penalty_hessian
%         Dim: problem size

import casadi.*

%% initialize NLP variable (stagewise, capital)
% define auxiliary variable dimension
eta_Dim = OCP.Dim.lambda;
% initialize problem variable 
X = SX.sym('X', OCP.Dim.x, OCP.nStages); 
XPrev = [OCP.x0, X(:, 1 : end - 1)];
U = SX.sym('U', OCP.Dim.u, OCP.nStages);
LAMBDA = SX.sym('LAMBDA', OCP.Dim.lambda, OCP.nStages);
ETA = SX.sym('ETA', eta_Dim, OCP.nStages); 
% initialize problem parameter
mu = SX.sym('mu', 1, 1);

%% mapping function object
% stage cost
L_S_map = OCP.FuncObj.L_S.map(OCP.nStages);
% D gap function
D_gap_map = self.D_gap_func.map(OCP.nStages);
% ODE r.h.s function
f_map = OCP.FuncObj.f.map(OCP.nStages);
% complementarity function
g_map = OCP.FuncObj.g.map(OCP.nStages);

%% formulate NLP function (stagewise)
% ocp stage cost
L_S_stage = L_S_map(X, U, LAMBDA);
% D gap
D_gap = D_gap_map(LAMBDA, ETA);
% constraint
f_stage = f_map(X, U, LAMBDA);
g_stage = g_map(X, U, LAMBDA);

%% reshape NLP variable and function (column, lowercase)
% variable
Z = [X; U; LAMBDA; ETA];
z = reshape(Z, (OCP.Dim.x + OCP.Dim.u + OCP.Dim.lambda + eta_Dim) * OCP.nStages, 1);
Dim.z_Node = cumsum([OCP.Dim.x, OCP.Dim.u, OCP.Dim.lambda, eta_Dim]); 
Dim.z = size(z, 1);
% problem parameter
p = mu;
% cost function
J_ocp = sum(L_S_stage) * OCP.timeStep;
J_penalty = mu * self.Huber_func(D_gap);
J = J_ocp + J_penalty;
% equality constraint h = 0
h_stage = ...
    [XPrev - X + f_stage * OCP.timeStep;...
    g_stage - ETA];
h = reshape(h_stage, (OCP.Dim.x + OCP.Dim.lambda) * OCP.nStages, 1);
Dim.h_Node = cumsum([OCP.Dim.x, OCP.Dim.lambda]);
Dim.h = size(h, 1);

%% formulate NLP jacobian and hessian
% cost gradient
J_grad = jacobian(J, z);

% constraint jacobian (constant matrix)
h_grad_formula = jacobian(h, z);
h_grad_func = Function('h_grad_func', {z}, {h_grad_formula}, {'z'}, {'h_grad_formula'});
h_grad = sparse(h_grad_func(zeros(Dim.z, 1)));

% ocp cost Hessian (constant matrix)
[J_ocp_hessian_formula, ~] = hessian(J_ocp, z);
J_ocp_hessian_func = Function('J_ocp_hessian_func', {z}, {J_ocp_hessian_formula}, {'z'}, {'J_ocp_hessian_formula'});
J_ocp_hessian = sparse(J_ocp_hessian_func(zeros(Dim.z, 1)));

D_gap_grad = jacobian(D_gap, z);

%% create output struct
nlp = struct('z', z, 'p', p,...
    'J', J, 'h', h, ...
    'J_ocp', J_ocp, 'J_penalty', J_penalty,...
    'D_gap', D_gap,...
    'J_grad', J_grad,...
    'h_grad', h_grad,...
    'J_ocp_hessian', J_ocp_hessian,... 
    'D_gap_grad', D_gap_grad,...
    'Dim', Dim);

end