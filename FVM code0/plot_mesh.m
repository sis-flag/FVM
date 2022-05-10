function plot_mesh(Mesh)

hold on

for E = 1:Mesh.nE
    np = Mesh.E2P{E};
    P0 = Mesh.xp{np(1)};
    P1 = Mesh.xp{np(2)};
    plot([P0(1), P1(1)], [P0(2), P1(2)], 'k-')
end

% for P = 1:Mesh.nP
%     txp = Mesh.xp{P};
%     text(txp(1), txp(2), num2str(P), 'Color', 'r')
% end
% for U = 1:Mesh.nU
%     txc = Mesh.xc{U};
%     text(txc(1), txc(2), num2str(U), 'Color', 'm')
% end
% for E = 1:Mesh.nE
%     txe = Mesh.xe{E};
%     text(txe(1), txe(2), num2str(E), 'Color', 'b')
% end
% for E = 1:Mesh.nE
%     quiver(Mesh.xe{E}(1), Mesh.xe{E}(2), ...
%         Mesh.nv{E}(1) / Mesh.nU, Mesh.nv{E}(2) / Mesh.nU, ...
%         'Color', 'k') 
% end