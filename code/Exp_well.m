clear;
addpath("mesh", "PDE", "scheme", "interp", "utils")

%% well 1
Mesh = load_mesh_file("mesh/mesh_well1.txt");

% figure
% plot_mesh(Mesh)
% set(gcf, 'Position', [200, 200, 330, 330])
% axis("equal")
% xlim([0, 1])
% ylim([0, 1])
% print('pics/w1m.eps','-depsc','-loose')

PDE = linear_stab_well1();

[A, F] = ECS1(Mesh, PDE);
u1 = A \ F;
fprintf("ECS-FV: max %g, min %g\n", max(u1), min(u1))

% figure
% temp_plot(Mesh, u1, 10, 12)
% set(gcf, 'Position', [200, 200, 400, 330])
% print('pics/w1s1.eps','-depsc','-loose')

[A, F] = ECS2(Mesh, PDE);
u2 = A \ F;
fprintf("ECS-MFD: max %g, min %g\n", max(u2), min(u2))

% figure
% temp_plot(Mesh, u2, 10, 12)
% set(gcf, 'Position', [200, 200, 400, 330])
% print('pics/w1s2.eps','-depsc','-loose')

%% well 2
Mesh = load_mesh_file("mesh/mesh_well2.txt");

% figure
% plot_mesh(Mesh)
% axis("equal")
% xlim([0, 1])
% ylim([0, 1])
% set(gcf, 'Position', [200, 200, 330, 330])
% print('pics/w2m.eps','-depsc','-loose')

PDE = linear_stab_well2();

[A, F] = ECS1(Mesh, PDE);
for E = 1:Mesh.nE
    xe = Mesh.xe(E); ye = Mesh.ye(E);
    if length(Mesh.E2U{E}) == 1 && ...
       0.1 < xe && xe < 0.9 && 0.1 < ye && ye < 0.9
        A(E, :) = 0;
        A(E, E) = 1;
        if xe < 0.5
            F(E) = 1;
        else
            F(E) = 0;
        end
    end
end
u1 = A \ F;
fprintf("ECS-FV: max %g, min %g\n", max(u1), min(u1))

% figure
% temp_plot(Mesh, u1, 0, 1)
% set(gcf, 'Position', [200, 200, 400, 330])
% print('pics/w2s1.eps','-depsc','-loose')

[A, F] = ECS2(Mesh, PDE);
for E = 1:Mesh.nE
    xe = Mesh.xe(E); ye = Mesh.ye(E);
    if length(Mesh.E2U{E}) == 1 && ...
       0.1 < xe && xe < 0.9 && 0.1 < ye && ye < 0.9
        A(E, :) = 0;
        A(E, E) = 1;
        if xe < 0.5
            F(E) = 1;
        else
            F(E) = 0;
        end
    end
end
u2 = A \ F;
fprintf("ECS-MFD: max %g, min %g\n", max(u2), min(u2))

% figure
% temp_plot(Mesh, u2, 0, 1)
% set(gcf, 'Position', [200, 200, 400, 330])
% print('pics/w2s2.eps','-depsc','-loose')

%% plot solutions, points out of range are marked
function temp_plot(Mesh, u, mu, Mu)
cmap = colormap();
ncolor = size(cmap, 1);
temp = (u - mu) / (Mu - mu) * (ncolor-1) + 1;
color = interp1(1:ncolor, cmap, temp, 'linear');

ind_out = (u < mu) | (u > Mu);
color(ind_out, :) = 0;

hold on
for E = 1:Mesh.nE
    nP = Mesh.E2P{E}; nU = Mesh.E2U{E};
    x1 = Mesh.xp(nP(1)); y1 = Mesh.yp(nP(1));
    x2 = Mesh.xp(nP(2)); y2 = Mesh.yp(nP(2));
    if length(nU) == 2
        x3 = Mesh.xc(nU(1)); y3 = Mesh.yc(nU(1));
        x4 = Mesh.xc(nU(2)); y4 = Mesh.yc(nU(2));
        xs = [x1, x3, x2, x4]; ys = [y1, y3, y2, y4];
    elseif length(nU) == 1
        x3 = Mesh.xc(nU); y3 = Mesh.yc(nU);
        xs = [x1, x3, x2]; ys = [y1, y3, y2];
    end

    if mu <= u(E) && u(E) <= Mu
        fill(xs, ys, color(E,:), 'LineStyle', 'none')
    else
        % fill(xs, ys, 'white', 'LineStyle', 'none')
    end

    plot([x1, x2], [y1, y2], 'k-', "LineWidth", 0.5)
end

colorbar
clim([mu, Mu])
axis("equal")
xlim([0, 1])
ylim([0, 1])
end
