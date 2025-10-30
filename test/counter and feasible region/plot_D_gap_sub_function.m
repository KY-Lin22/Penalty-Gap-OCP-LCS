clear all
clc

counter_level = [-10, -1, 1, 10];

plotFBSQRTFunContour(counter_level)

%%
function plotFBSQRTFunContour(counter_level)
% parameter
a = 0.5;
b = 2;
stepsize = 0.01;
% data
x = -10 : stepsize : 10;
y = -10 : stepsize : 10;
Z_surf = zeros(length(x), length(y));
for i = 1 : length(x)
    for j = 1 : length(y)
        Z_surf(j, i) = -a/2 * x(i)^2  + x(i) * y(j) - 1/(2*b) * y(j)^2;
        % Z_surf(j, i) = b/2 * x(i)^2  - x(i) * y(j) + 1/(2*a) * y(j)^2;
    end
end
% mesh and colour
[X, Y] = meshgrid(x, y);
C = 1.*Z_surf;
figure(1)
sfc_handle = surfc(X, Y, Z_surf, C, 'FaceAlpha',0.5, 'EdgeColor', 'none');

contourProperty = sfc_handle(2);
contourProperty.ContourZLevel = -20;
contourProperty.LineWidth = 1.5;
contourProperty.LevelList = counter_level;
contourProperty.ShowText = 'on';
hold on

colorbar
% contour
% surf
% fsurf
lighting gouraud;% best lighting algoithm for curved surfaces
material shiny
box on
% sets the axis limits equal to the range of the data
axis([x(1)-0.5, x(end)+0.5, y(1)-0.5, y(end)+0.5, contourProperty.ContourZLevel, 25] )
set(gca, 'FontSize', 20)
xlabel('$\lambda_i$', 'Interpreter','latex', 'FontSize', 20)
ylabel('$\eta_i$', 'Interpreter','latex', 'FontSize', 20)
zlabel('$\delta^{ab}$', 'Interpreter','latex', 'FontSize', 20)
end