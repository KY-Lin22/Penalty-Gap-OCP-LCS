clear all
clc

a = 0.7;
b = 1.5;

plot_D_gap_function_convex_region(a, b);

%%
function plot_D_gap_function_convex_region(a, b)

% parameter
stepsize = 0.01;

% figure limit
x_lb = -8;
x_ub = 8;
y_lb = -8;
y_ub = 8;

% node (x, y)
origin = [0, 0];
node_1 = [0, y_ub];
node_2 = [y_ub/b, y_ub];
node_3 = [x_ub, y_ub];
node_4 = [x_ub, a*x_ub];
node_5 = [x_ub, 0];
node_6 = [x_ub, y_lb];
node_7 = [y_lb/b, y_lb];
node_8 = [x_lb, y_lb];
node_9 = [x_lb, a*x_lb];
node_10 = [x_lb, y_ub];

% axis
bound_y_ax_x = x_lb: stepsize : x_ub;
bound_y_ax_y = a.* bound_y_ax_x;
bound_y_bx_x = x_lb: stepsize : x_ub;
bound_y_bx_y = b.* bound_y_bx_x;

% plot boundary
figure(2)
plot(bound_y_ax_x, bound_y_ax_y, 'k', 'LineWidth', 2)
hold on
plot(bound_y_bx_x, bound_y_bx_y, 'k', 'LineWidth', 2)
hold on
plot([0, 0], [0, y_ub], 'k', 'LineWidth', 2)
hold on
plot([0, x_ub], [0, 0], 'k', 'LineWidth', 2)
hold on
% color region
patch([origin(1), node_1(1), node_2(1)], [origin(2), node_1(2), node_2(2)],...
    'blue', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
patch([origin(1), node_2(1), node_3(1), node_4(1)], [origin(2), node_2(2), node_3(2), node_4(2)],...
    'green', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
patch([origin(1), node_4(1), node_5(1)], [origin(2), node_4(2), node_5(2)],...
    'yellow', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
patch([origin(1), node_5(1), node_6(1), node_7(1)], [origin(2), node_5(2), node_6(2), node_7(2)],...
    'red', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
patch([origin(1), node_7(1), node_8(1), node_9(1)], [origin(2), node_7(2), node_8(2), node_9(2)],...
    'cyan', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
patch([origin(1), node_9(1), node_10(1), node_1(1)], [origin(2), node_9(2), node_10(2), node_1(2)],...
    'magenta', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
% text axis
text(1, -1, '$\eta_i = a \lambda_i$', 'Interpreter','latex', 'FontSize', 20)
annotation('arrow', [0.1, 0.3], [0.2, 0.4], 'HeadWidth', 10, 'HeadLength', 10, 'LineWidth', 1.2)
text(-4, 1, '$\eta_i = b \lambda_i$', 'Interpreter','latex', 'FontSize', 20)
annotation('arrow', [0.4, 0.6], [0.7, 0.8], 'HeadWidth', 10, 'HeadLength', 10, 'LineWidth', 1.2)
% text multiplier
text(0.1, 7, '$v_i < 0, w_i = 0$', 'Interpreter','latex', 'FontSize', 15)
text(4.1, 5.7, '$v_i < 0,$', 'Interpreter','latex', 'FontSize', 15)
text(4.1, 4.7, '$w_i < 0$', 'Interpreter','latex', 'FontSize', 15)
text(3, 1, '$v_i = 0, w_i < 0$', 'Interpreter','latex', 'FontSize', 15)
text(2, -4, '$v_i = 0, w_i > 0$', 'Interpreter','latex', 'FontSize', 15)
text(-6.9, -5.8, '$v_i > 0,$', 'Interpreter','latex', 'FontSize', 15)
text(-6.9, -6.8, '$w_i > 0$', 'Interpreter','latex', 'FontSize', 15)
text(-6, 4, '$v_i > 0, w_i = 0$', 'Interpreter','latex', 'FontSize', 15)

axis([x_lb, x_ub, y_lb, y_ub])
set(gca, 'FontSize', 20)
xline(0, 'LineWidth', 1); 
yline(0, 'LineWidth', 1);
xlabel('$\lambda_i$', 'Interpreter','latex', 'FontSize', 20)
ylabel('$\eta_i$', 'Interpreter','latex', 'FontSize', 20)



end