function plot_mesh(Mesh)

hold on

for E = 1:Mesh.nE
    nPs = Mesh.E2P{E};
    P1 = nPs(1); P2 = nPs(2);
    x1 = Mesh.xp(P1); x2 = Mesh.xp(P2);
    y1 = Mesh.yp(P1); y2 = Mesh.yp(P2);
    plot([x1, x2], [y1, y2], 'k-')
end

% mark the numbers of each element
if Mesh.nP < 40
    for P = 1:Mesh.nP
        xp = Mesh.xp(P); yp = Mesh.yp(P);
        text(xp, yp, num2str(P), 'Color', 'r')
    end
    for U = 1:Mesh.nU
        xc = Mesh.xc(U); yc = Mesh.yc(U);
        text(xc, yc, num2str(U), 'Color', 'm')
    end
    for E = 1:Mesh.nE
        xe = Mesh.xe(E); ye = Mesh.ye(E);
        text(xe, ye, num2str(E), 'Color', 'b')

        nx = Mesh.nx(E) / Mesh.nU;
        ny = Mesh.ny(E) / Mesh.nU;
        quiver(xe, ye, nx, ny, 'Color', 'k')
    end
end