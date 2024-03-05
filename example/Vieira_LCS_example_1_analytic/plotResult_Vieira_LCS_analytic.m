function plotResult_Vieira_LCS_analytic(OCP, NLP, z_Opt)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

Z_Opt = reshape(z_Opt, NLP.Dim.z_Node(end), OCP.nStages);

X_Opt = Z_Opt(1 : NLP.Dim.z_Node(1), :);
U_Opt = Z_Opt(NLP.Dim.z_Node(1) + 1 : NLP.Dim.z_Node(2), :);
LAMBDA_Opt = Z_Opt(NLP.Dim.z_Node(2) + 1 : NLP.Dim.z_Node(3), :);

timeAxis = 0 : OCP.timeStep : OCP.nStages * OCP.timeStep;

g_FuncObj_map = OCP.FuncObj.g.map(OCP.nStages);
g_Opt = full(g_FuncObj_map(X_Opt, U_Opt, LAMBDA_Opt));

%% state, control, and algebraic trajectory
figure(111)
subplot(2,2,1)
plot(timeAxis, [OCP.x0, X_Opt], 'r', 'LineWidth',1.2)
xlabel('time [s]')
title('state')

subplot(2,2,2)
plot(timeAxis(2:end), LAMBDA_Opt, 'r', 'LineWidth', 1.2, 'LineStyle', 'none', 'Marker', '.')
xlabel('time [s]')
title('algebraic var')

subplot(2,2,3)
stairs(timeAxis(2:end), U_Opt, 'LineWidth', 1.2)
xlabel('time [s]')
title('control')

subplot(2,2,4)
plot(timeAxis(2:end), g_Opt, 'b', 'LineWidth', 1.2, 'LineStyle', 'none', 'Marker', '.')
xlabel('time [s]')
title('VI func')

%% anlytical optimal trajectory 
[X_analytic_Opt, U_analytic_Opt, LAMBDA_analytic_Opt] = compute_analytic_optimal_trajectory(OCP);
g_analytic_Opt = full(g_FuncObj_map(X_analytic_Opt, U_analytic_Opt, LAMBDA_analytic_Opt));

figure(112)
subplot(2,2,1)
plot(timeAxis, [OCP.x0, X_analytic_Opt], 'r', 'LineWidth',1.2)
xlabel('time [s]')
title('state (analytical)')
subplot(2,2,2)
plot(timeAxis(2:end), LAMBDA_analytic_Opt, 'r', 'LineWidth', 1.2, 'LineStyle', 'none', 'Marker', '.')
xlabel('time [s]')
title('algebraic var (analytical)')
subplot(2,2,3)
stairs(timeAxis(2:end), U_analytic_Opt, 'LineWidth', 1.2)
xlabel('time [s]')
title('control (analytical)')
subplot(2,2,4)
plot(timeAxis(2:end), g_analytic_Opt, 'b', 'LineWidth', 1.2, 'LineStyle', 'none', 'Marker', '.')
xlabel('time [s]')
title('VI func (analytical)')

%% trajectory error

figure(113)
subplot(2,2,1)
plot(timeAxis(2:end), abs(X_Opt - X_analytic_Opt), 'r', 'LineWidth',1.2)
xlabel('time [s]')
title('state error')

subplot(2,2,2)
plot(timeAxis(2:end), abs(LAMBDA_Opt - LAMBDA_analytic_Opt), 'r', 'LineWidth', 1.2, 'LineStyle', 'none', 'Marker', '.')
xlabel('time [s]')
title('algebraic var error')

subplot(2,2,3)
stairs(timeAxis(2:end), abs(U_Opt - U_analytic_Opt), 'LineWidth', 1.2)
xlabel('time [s]')
title('control error')

subplot(2,2,4)
plot(timeAxis(2:end), abs(g_Opt - g_analytic_Opt), 'b', 'LineWidth', 1.2, 'LineStyle', 'none', 'Marker', '.')
xlabel('time [s]')
title('VI func error')


end