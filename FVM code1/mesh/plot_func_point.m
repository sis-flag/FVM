function plot_func_point(Mesh, up)

hold on

cmap = colormap('jet');
temp = (up - min(up)) / range(up) * (size(cmap,1)-1);
color = cmap(ceil(temp + 0.5), :);

for E = 1:Mesh.nE
    nPs = Mesh.E2P{E};
    P1 = nPs(1); P2 = nPs(2);
    x1 = Mesh.xp(P1); x2 = Mesh.xp(P2);
    y1 = Mesh.yp(P1); y2 = Mesh.yp(P2);
    u1 = up(P1); u2 = up(P2);
    plot3([x1, x2], [y1, y2], [u1, u2], ...
        'Color', (color(P1,:) + color(P2,:))/2)
end


end