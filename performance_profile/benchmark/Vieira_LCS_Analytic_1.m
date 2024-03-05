function OCP = Vieira_LCS_Analytic_1()
% ref: example 1 in ''Quadratic Optimal Control of Linear Complementarity 
%      Systems : First order necessary conditions and numerical analysis''
%      2018, A. Vieira, et. al,
import casadi.*

% time parameter
TimeHorizon = 1; % time horizon T
nStages = 100; % number of discretized stages
timeStep = TimeHorizon ./ nStages; % discretization time step

% initial state
x0 = 1; 

% variable
xDim = 1;
uDim = 1;
lambdaDim = 1;
x = SX.sym('x', xDim, 1);
u = SX.sym('u', uDim, 1);
lambda = SX.sym('lambda', lambdaDim, 1);

% cost function
Q_x = 1;
Q_u = 1;
Q_lambda = 0;

% LCS
param.a = 3;
param.b = -0.5;
param.d = 1;
param.e = -2;
param.f = 3;

A = param.a;
B = param.f;
E = param.b;

C = 0;
D = param.e;
F = param.d;

%% create OCP instant
OCP = OCP_Formulation(...
    TimeHorizon, nStages, timeStep,...
    x0, ...
    x, u, lambda,...
    Q_x, Q_u, Q_lambda,...
    A, B, E,...
    C, D, F);
end