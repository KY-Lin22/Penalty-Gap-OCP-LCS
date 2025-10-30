function OCP = OCP_Vieira_LCS_state_jump()
% ref: example 7 in ''Quadratic Optimal Control of Linear Complementarity 
%      Systems : First order necessary conditions and numerical analysis''
%      2018, A. Vieira, et. al,
import casadi.*

% time parameter
TimeHorizon = 10; % time horizon T
nStages = 500; % number of discretized stages
timeStep = TimeHorizon ./ nStages; % discretization time step

x0 = [-2; 1; -1]; % initial state

% variable
xDim = 3;
uDim = 1;
lambdaDim = 1;
x = SX.sym('x', xDim, 1);
u = SX.sym('u', uDim, 1);
lambda = SX.sym('lambda', lambdaDim, 1);

% cost function
lambdaWeight = 0.1;
Q_x = diag([1, 1, 1]);
Q_u = 1;
Q_lambda = lambdaWeight;

% LCS
A = [0, 1, 0; 0, 0, 1; 0, 0, 0];
B = [0; 0; 1];
E = [0; 0; 1];

C = [1, 0, 0];
D = 1;
F = 0;

%% create OCP instant
OCP = OCP_Formulation(...
    TimeHorizon, nStages, timeStep,...
    x0, ...
    x, u, lambda,...
    Q_x, Q_u, Q_lambda,...
    A, B, E,...
    C, D, F);

end