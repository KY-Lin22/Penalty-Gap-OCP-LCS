function OCP = Vieira_LCS_With_Penalty_1()
% ref: example 5 in ''Quadratic Optimal Control of Linear Complementarity 
%      Systems : First order necessary conditions and numerical analysis''
%      2018, A. Vieira, et. al,
import casadi.*

% time parameter
TimeHorizon = 1; % time horizon T
nStages = 100; % number of discretized stages
timeStep = TimeHorizon ./ nStages; % discretization time step

x0 = [-0.5; -1]; % initial state

% variable
xDim = 2;
uDim = 1;
lambdaDim = 1;
x = SX.sym('x', xDim, 1);
u = SX.sym('u', uDim, 1);
lambda = SX.sym('lambda', lambdaDim, 1);

% cost function
alpha = 10;
Q_x = diag([1, 1]);
Q_u = 1;
Q_lambda = alpha;

% LCS
A = [5, -6; 3, 9];
B = [0; -4];
E = [4; 5];

C = [-1, 5];
D = 6;
F = 1;

%% create OCP instant
OCP = OCP_Formulation(...
    TimeHorizon, nStages, timeStep,...
    x0, ...
    x, u, lambda,...
    Q_x, Q_u, Q_lambda,...
    A, B, E,...
    C, D, F);
end