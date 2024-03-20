%% counter of D gap function
clear all
clc

a = 0.5;
b = 2;
counter_level = [1, 5, 10];
plot_D_gap_function_convex_region(a, b, counter_level);

%%
function plot_D_gap_function_convex_region(a, b, counter_level)
% parameter
stepsize = 0.01;
num_level = numel(counter_level);

% node
origin_x = 0;
origin_y = 0;

nodePoint1_x = zeros(1, num_level);
nodePoint2_x = zeros(1, num_level);
nodePoint2_y = zeros(1, num_level);
nodePoint3_x = zeros(1, num_level);
nodePoint4_x = zeros(1, num_level);
nodePoint4_y = zeros(1, num_level);
for i = 1 : num_level
    s = counter_level(i);
    nodePoint1_x(1, i) = sqrt(2*s/(b-a));
    nodePoint2_x(1, i) = sqrt(2*b*s/(a*(b-a)));
    nodePoint2_y(1, i) = sqrt(2*a*b*s/(b-a));
    nodePoint3_x(1, i) = -sqrt(2*s/(b-a));
    nodePoint4_x(1, i) = -sqrt(2*a*s/(b*(b-a)));
    nodePoint4_y(1, i) = -sqrt(2*a*b*s/(b-a));
end
nodePoint1_y = b*nodePoint1_x;
nodePoint3_y = a*nodePoint3_x;

% figure boundary
x_lb = nodePoint3_x(end) - 1;
x_ub = nodePoint2_x(end) + 1;
y_lb = nodePoint4_y(end) - 1;
y_ub = nodePoint1_y(end) + 1;

% axis
bound_y_ax_x = x_lb: stepsize : x_ub;
bound_y_ax_y = a.* bound_y_ax_x;
bound_y_bx_x = x_lb: stepsize : x_ub;
bound_y_bx_y = b.* bound_y_bx_x;

% convex and concave region
convex_reg_up_x = [origin_x, bound_y_bx_x(end), x_lb, bound_y_ax_x(1)];
convex_reg_up_y = [origin_y, bound_y_bx_y(end), y_ub, bound_y_ax_y(1)];

convex_reg_down_x = [bound_y_ax_x(1), bound_y_ax_x(end), x_ub, x_lb];
convex_reg_down_y = [bound_y_ax_y(1), bound_y_ax_y(end), y_lb, y_lb];

concave_reg_x = [origin_x, bound_y_bx_x(end), x_ub, bound_y_ax_x(end)];
concave_reg_y = [origin_y, bound_y_bx_y(end), y_ub, bound_y_ax_y(end)];

% plot counter and boundary
figure(2)
for i = 1 : num_level
    s = counter_level(i);
    plot([nodePoint3_x(i), nodePoint3_x(i)], [y_ub, nodePoint3_y(i)], 'k', 'LineStyle', '--', 'LineWidth', 1)
    hold on
    plot([nodePoint1_x(i), nodePoint1_x(i)], [y_ub, nodePoint1_y(i)], 'k', 'LineStyle', '--', 'LineWidth', 1)
    hold on
    plot([nodePoint2_x(i), x_ub], [nodePoint2_y(i), nodePoint2_y(i)], 'k', 'LineStyle', '--', 'LineWidth', 1)
    hold on
    plot([nodePoint4_x(i), x_ub], [nodePoint4_y(i), nodePoint4_y(i)], 'k', 'LineStyle', '--', 'LineWidth', 1)
    hold on
    f2 = @(x,y) -a./2 .* x.^2 + x .* y - 1./(2.* b) .* y.^2 - s;
    fimplicit(f2, [nodePoint1_x(i), nodePoint2_x(i), nodePoint2_y(i), nodePoint1_y(i)], 'k', 'LineStyle', '--', 'LineWidth', 1) % region 2 curve
    hold on
    f4 = @(x,y) b./2 .* x.^2 - x .* y + 1./(2.* a) .* y.^2 - s;
    fimplicit(f4, [nodePoint3_x(i), nodePoint4_x(i), nodePoint4_y(i), nodePoint3_y(i)], 'k', 'LineStyle', '--', 'LineWidth', 1) % region 4 curve
    hold on
end
plot(bound_y_ax_x, bound_y_ax_y, 'k', 'LineWidth', 2)
hold on
plot(bound_y_bx_x, bound_y_bx_y, 'k', 'LineWidth', 2)
hold on

% plot concave region
patch(convex_reg_up_x, convex_reg_up_y, 'green', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
patch(convex_reg_down_x, convex_reg_down_y, 'green', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
patch(concave_reg_x, concave_reg_y, 'blue', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on


axis([nodePoint3_x(end)-0.6, nodePoint2_x(end)+0.5, nodePoint4_y(end)-0.4, nodePoint1_y(end)+0.5])
xline(0, 'LineWidth', 1); 
yline(0, 'LineWidth', 1);

end