function OCP = Vieira_LCS_High_Dim()
% ref: example 3 in ''Quadratic Optimal Control of Linear Complementarity 
%      Systems : First order necessary conditions and numerical analysis''
%      2018, A. Vieira, et. al,
import casadi.*

% time parameter
TimeHorizon = 1; % time horizon T
nStages = 100; % number of discretized stages
timeStep = TimeHorizon ./ nStages; % discretization time step

x0 = [-0.5; 1]; % initial state

% variable
xDim = 2;
uDim = 2;
lambdaDim = 2;
x = SX.sym('x', xDim, 1);
u = SX.sym('u', uDim, 1);
lambda = SX.sym('lambda', lambdaDim, 1);

% cost function
Q_x = diag([1, 1]);
Q_u = diag([25, 25]);
Q_lambda = diag([0, 0]);

A = [1, 2; 2, 1];
B = [1, 3; 2, 1];
E = [-1, 1; -1, 1];

C = [3, -1; -2, 0];
D = [1, -1; -1, 2];
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