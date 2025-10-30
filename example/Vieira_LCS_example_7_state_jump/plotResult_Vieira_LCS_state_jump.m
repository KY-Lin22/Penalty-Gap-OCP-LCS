function plotResult_Vieira_LCS_state_jump(OCP, NLP, z_Opt)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
Z_Opt = reshape(z_Opt, NLP.Dim.z_Node(end), OCP.nStages);

X_Opt = Z_Opt(1 : NLP.Dim.z_Node(1), :);
U_Opt = Z_Opt(NLP.Dim.z_Node(1) + 1 : NLP.Dim.z_Node(2), :);
LAMBDA_Opt = Z_Opt(NLP.Dim.z_Node(2) + 1 : NLP.Dim.z_Node(3), :);

timeAxis = 0 : OCP.timeStep : OCP.nStages * OCP.timeStep;

f_FuncObj_map = OCP.FuncObj.f.map(OCP.nStages);
f_Opt = full(f_FuncObj_map(X_Opt, U_Opt, LAMBDA_Opt));

g_FuncObj_map = OCP.FuncObj.g.map(OCP.nStages);
g_Opt = full(g_FuncObj_map(X_Opt, U_Opt, LAMBDA_Opt));

figure(111)
subplot(3,2,3)
plot(timeAxis, [OCP.x0(1), X_Opt(1, :)], 'r', 'LineWidth',1.2)
xlabel('time [s]', 'Interpreter','latex')
title('state $x_1$', 'Interpreter','latex')

subplot(3,2,4)
plot(timeAxis, [OCP.x0(2), X_Opt(2, :)], 'g', 'LineWidth',1.2)
xlabel('time [s]', 'Interpreter','latex')
title('state $x_2$', 'Interpreter','latex')

subplot(3,2,5)
plot(timeAxis(2:end), f_Opt(1, :), 'r-', 'LineWidth', 1.2)
xlabel('time [s]', 'Interpreter','latex')
title('derivative $\dot{x}_1$', 'Interpreter','latex')

subplot(3,2,6)
plot(timeAxis(2:end), f_Opt(2, :), 'g-', 'LineWidth', 1.2)
xlabel('time [s]', 'Interpreter','latex')
title('derivative $\dot{x}_2$', 'Interpreter','latex')

figure(112)
subplot(3,2,1)
plot(timeAxis, [OCP.x0(3), X_Opt(3, :)], 'b', 'LineWidth',1.2)
xlabel('time [s]', 'Interpreter','latex')
title('state $x_3$', 'Interpreter','latex')

subplot(3,2,2)
plot(timeAxis(2:end), LAMBDA_Opt(1, :), 'r', 'LineWidth', 1.2)
xlabel('time [s]', 'Interpreter','latex')
title('algebraic var. $ \lambda$', 'Interpreter','latex')

subplot(3,2,3)
plot(timeAxis(2:end), f_Opt(3, :), 'b-', 'LineWidth', 1.2)
xlabel('time [s]', 'Interpreter','latex')
title('derivative $\dot{x}_3$', 'Interpreter','latex')

subplot(3,2,4)
plot(timeAxis(2:end), g_Opt(1, :), 'm', 'LineWidth', 1.2)
xlabel('time [s]', 'Interpreter','latex')
title('function $g$', 'Interpreter','latex')

subplot(3,2,5)
stairs(timeAxis(2:end), U_Opt(1,:), 'LineWidth', 1.2)
xlabel('time [s]', 'Interpreter','latex')
title('control $u$', 'Interpreter','latex')

end

