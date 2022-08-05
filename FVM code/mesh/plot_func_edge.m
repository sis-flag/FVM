function plot_func_edge(Mesh, ue)

hold on

cmap = colormap('jet');
temp = (ue - min(ue)) / range(ue) * (size(cmap,1)-1);
color = cmap(ceil(temp + 0.5), :);

for E = 1:Mesh.nE
    nP = Mesh.E2P{E};
    x1 = Mesh.xp(nP(1)); y1 = Mesh.yp(nP(1));
    x2 = Mesh.xp(nP(2)); y2 = Mesh.yp(nP(2));
    
    u0 = ue(E);
    plot3([x1, x2], [y1, y2], [u0, u0], ...
        'Color', color(E,:), 'LineWidth', 1)
end

end