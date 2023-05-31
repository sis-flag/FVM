function pcolor_func_edge(Mesh, ue)

hold on

Mu = max(ue);
mu = min(ue);

cmap = colormap();
ncolor = size(cmap, 1);
temp = (ue - mu) / (Mu - mu) * (ncolor-1) + 1;
color = interp1(1:ncolor, cmap, temp, 'linear');

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
    
    fill(xs, ys, color(E,:), 'LineStyle', 'none')
    plot([x1, x2], [y1, y2], 'k-')
end

colorbar
clim([mu, Mu])
end