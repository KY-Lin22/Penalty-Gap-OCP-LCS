function natRes = evaluate_natural_residual(self, z_Opt)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*

NLP = self.NLP;
OCP = self.OCP;

%% extract solution
Z_Opt = reshape(z_Opt, NLP.Dim.z_Node(end), OCP.nStages);

X_Opt = Z_Opt(1 : NLP.Dim.z_Node(1), :);
U_Opt = Z_Opt(NLP.Dim.z_Node(1) + 1 : NLP.Dim.z_Node(2), :);
if strcmp(NLP.penalty_problem, 'gap_polar_based')
    THETA_Opt = Z_Opt(NLP.Dim.z_Node(2) + 1 : NLP.Dim.z_Node(3), :);
    R_Opt = Z_Opt(NLP.Dim.z_Node(3) + 1 : NLP.Dim.z_Node(4), :);
    % polar function
    polar_func_map = NLP.polar_func.map(OCP.nStages);
    [LAMBDA_Opt, ~] = polar_func_map(THETA_Opt, R_Opt);
    LAMBDA_Opt = full(LAMBDA_Opt);
else
    LAMBDA_Opt = Z_Opt(NLP.Dim.z_Node(2) + 1 : NLP.Dim.z_Node(3), :);
end


%% problem data for Euclidean projector
g_FuncObj_map = OCP.FuncObj.g.map(OCP.nStages);
g_Opt = full(g_FuncObj_map(X_Opt, U_Opt, LAMBDA_Opt));
w_Opt = LAMBDA_Opt - g_Opt; 

%% construct qp solver to evaluate Euclidean projector
% variable(projector output) and parameter(projector input)
y = SX.sym('y', OCP.Dim.lambda, 1);
w = SX.sym('w', OCP.Dim.lambda, 1);
% cost function and constraint (y >= 0)
J = 0.5 * y'*diag(ones(OCP.Dim.lambda, 1))*y - w'*y;
lbg = zeros(OCP.Dim.lambda, 1);
ubg = inf*ones(OCP.Dim.lambda, 1);
% problem struct
Prob = struct('x', y, 'f', J, 'g', y, 'p', w);
% option
Option = struct(...
    'printLevel', 'none',...% 'none', 'low', 'medium', 'high' (see qpoases manual sec 5.2)
    'hessian_type', 'posdef',...% 'unknown', 'posdef', 'semidef', 'indef', 'zero', 'identity' (see qpoases manual sec 4.4, 4.5)
    'error_on_fail', false);
solver_singleStage = qpsol('EuclideanProjector', 'qpoases', Prob, Option);
% solver
EuclideanProjector = solver_singleStage.map(OCP.nStages);

solution = EuclideanProjector(...
    'lbg', repmat(lbg, 1, OCP.nStages),...
    'ubg', repmat(ubg, 1, OCP.nStages),...
    'p', w_Opt);
proj_w_Opt = full(solution.x);

%% evaluate natural residual
% natRes = sum(sum(abs(LAMBDA_Opt - proj_w_Opt))) / (OCPEC.nStages * OCPEC.Dim.lambda); % sum of all natRes and then scaling
natRes = max(max(abs(LAMBDA_Opt - proj_w_Opt)));
end