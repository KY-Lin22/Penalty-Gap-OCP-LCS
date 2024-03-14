%% counter of D gap function
clear all
clc

a = 0.5;
b = 2;
plotDGapFunContour(a, b)

%%
function plotDGapFunContour(a, b)
% parameter
stepsize = 0.01;
% data
x = -4 : stepsize : 6;
y = -4 : stepsize : 6;
Z_surf = zeros(length(x), length(y));
for i = 1 : length(x)
    for j = 1 : length(y)
        Z_surf(j, i) = (b-a)/(2*a*b)*y(j)^2 - 1/(2*a)*(max([0, y(j)-a*x(i)]))^2 + 1/(2*b)*(max([0,y(j)-b*x(i)]))^2;
    end
end
% mesh and colour
[X, Y] = meshgrid(x, y);
C = 1.*Z_surf;
figure(3)
sfc_handle = surfc(X, Y, Z_surf, C, 'FaceAlpha',0.5, 'EdgeColor', 'none');
contourProperty = sfc_handle(2);
contourProperty.ContourZLevel = -max(max(Z_surf));
contourProperty.LineWidth = 1;
contourProperty.LevelList = [1, 5, 10];
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
axis([x(1)-0.5, x(end)+0.5, y(1)-0.5, y(end)+0.5, contourProperty.ContourZLevel, max(max(Z_surf)) + 5] )
% title('generalized D Gap Function')
xlabel('$\lambda_i$', 'Interpreter','latex')
ylabel('$\eta_i$', 'Interpreter','latex')
zlabel('$\delta^{ab}$', 'Interpreter','latex')
end