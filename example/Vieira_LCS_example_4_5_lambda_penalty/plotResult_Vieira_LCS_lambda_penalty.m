function plotResult_Vieira_LCS_lambda_penalty(OCP, NLP, z_Opt)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
Z_Opt = reshape(z_Opt, NLP.Dim.z_Node(end), OCP.nStages);

X_Opt = Z_Opt(1 : NLP.Dim.z_Node(1), :);
U_Opt = Z_Opt(NLP.Dim.z_Node(1) + 1 : NLP.Dim.z_Node(2), :);
LAMBDA_Opt = Z_Opt(NLP.Dim.z_Node(2) + 1 : NLP.Dim.z_Node(3), :);

timeAxis = 0 : OCP.timeStep : OCP.nStages * OCP.timeStep;

g_FuncObj_map = OCP.FuncObj.g.map(OCP.nStages);
g_Opt = full(g_FuncObj_map(X_Opt, U_Opt, LAMBDA_Opt));

figure(111)
subplot(2,2,1)
plot(timeAxis, [OCP.x0(1), X_Opt(1, :)], 'r',...
     timeAxis, [OCP.x0(2), X_Opt(2, :)], 'g', 'LineWidth',1.2)
legend('x1', 'x2')
xlabel('time [s]')
title('differential state')

subplot(2,2,2)
plot(timeAxis(2:end), LAMBDA_Opt(1, :), 'r', 'LineWidth', 1.2, 'LineStyle', 'none', 'Marker', '.')
xlabel('time [s]')
title('algebraic variable')

subplot(2,2,3)
stairs(timeAxis(2:end), U_Opt(1,:), 'LineWidth', 1.2)
xlabel('time [s]')
title('control input')

subplot(2,2,4)
plot(timeAxis(2:end), g_Opt(1, :), 'b', 'LineWidth', 1.2, 'LineStyle', 'none', 'Marker', '.')
xlabel('time [s]')
title('VI function')

end

