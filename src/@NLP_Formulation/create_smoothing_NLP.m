function nlp = create_smoothing_NLP(self, OCP)
% formulate a smoothing based NLP
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
%       c(z, p) >= 0
% where: (1) z: collects all the variables to be optimized and arranged in a stagewise manner
%            z = [z_1;...z_n;...z_N] and z_n = [x_n; u_n; lambda_n; eta_n] 
%            with x_n:      system state
%                 u_n:      control                 
%                 lambda_n: algebraic variable   
%                 eta_n:  auxiliary variable for complementarity function
%        (2) p: collects all the NLP problem parameters p = [mu]
%            with mu: nonnegative smoothing parameter
%        (3) J: cost function J = J_ocp + J_penalty
%            J_ocp = sum(J_stage_n)*dt
%            J_penalty = 0
%        (4) h: equality constraint arranged in a stagewise manner
%            h = [h_1;...h_n;...h_N] and 
%            h_n = [x_{n-1} - x_n + f(x_n, u_n, lambda_n)*dt;
%                   g(x_n, u_n, lambda_n) - eta_n;...
%                   FB(lambda_n, eta_n, 1/mu)]     
%        (5) c: this reformulation does not have inequality constraint
%
% output: nlp is a structure with fields:
%         z, p: variable and parameter
%         J, h, c: cost and constraint function,  
%         J_ocp, J_penalty: cost term
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

%% smoothing function
switch self.smoothing_problem
    case 'FB'
        u = SX.sym('u', OCP.Dim.lambda, 1);
        v = SX.sym('v', OCP.Dim.lambda, 1);
        s = SX.sym('s', 1, 1); % smoothing parameter
        FB_func = sqrt(u.^2 + v.^2 + s^2) - u - v;
        smoothing_func = Function('smoothing_func', {u, v, s}, {FB_func}, {'u', 'v', 's'}, {'FB_func'});
end

%% mapping function object
% stage cost
L_S_map = OCP.FuncObj.L_S.map(OCP.nStages);
% ODE r.h.s function
f_map = OCP.FuncObj.f.map(OCP.nStages);
% complementarity function
g_map = OCP.FuncObj.g.map(OCP.nStages);
% smoothing function
smoothing_func_map = smoothing_func.map(OCP.nStages);

%% formulate NLP function (stagewise)
% stage cost
L_S_stage = L_S_map(X, U, LAMBDA);
% constraint
f_stage = f_map(X, U, LAMBDA);
g_stage = g_map(X, U, LAMBDA);
smooth_func_stage = smoothing_func_map(LAMBDA, ETA, 1/mu);

%% reshape NLP variable and function (column, lowercase)
% variable
Z = [X; U; LAMBDA; ETA];
z = reshape(Z, (OCP.Dim.x + OCP.Dim.u + OCP.Dim.lambda + eta_Dim) * OCP.nStages, 1);
Dim.z_Node = cumsum([OCP.Dim.x, OCP.Dim.u, OCP.Dim.lambda, eta_Dim]); 
Dim.z = size(z, 1);
% problem parameter
p = mu;
Dim.p = size(p, 1);

% cost function
J_ocp = sum(L_S_stage)*OCP.timeStep;
J_penalty = 0;
J = J_ocp + J_penalty;

% equality constraint h = 0
h_stage = ...
    [XPrev - X + f_stage * OCP.timeStep;...
    g_stage - ETA;...
    smooth_func_stage];
h = reshape(h_stage, (OCP.Dim.x + OCP.Dim.lambda + OCP.Dim.lambda) * OCP.nStages, 1);
Dim.h_Node = cumsum([OCP.Dim.x, OCP.Dim.lambda, OCP.Dim.lambda]);
Dim.h = size(h, 1);

% inequality constraint c >= 0
c = SX.sym('c', 0, 1);
Dim.c = size(c, 1);

%% create output struct
nlp = struct('z', z, 'p', p,...
    'J', J, 'h', h, 'c', c, ... 
    'J_ocp', J_ocp, 'J_penalty', J_penalty,...
    'Dim', Dim);

end