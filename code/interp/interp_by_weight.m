function up = interp_by_weight(Mesh, PDE, uc, weight)
% get interplation weights on mesh points by average interplation

up = zeros(Mesh.nP, 1);
for P = 1:Mesh.nP
    if Mesh.isbdp(P)
        xp = Mesh.xp(P); yp = Mesh.yp(P);
        up(P) = PDE.u(xp, yp);
    else
        nU = Mesh.P2U{P};
        up(P) = weight{P}' * uc(nU);
    end
end

end