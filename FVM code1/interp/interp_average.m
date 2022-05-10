function up = interp_average(Mesh, PDE, uc)
% get interplation weights on mesh points by average interplation

up = zeros(Mesh.nP, 1);
for P = 1:Mesh.nP
    if Mesh.isbdp(P)
        xp = Mesh.xp(P); yp = Mesh.yp(P);
        up(P) = PDE.u(xp, yp);
    else
        nU = Mesh.P2U{P};
        up(P) = mean(uc(nU));
    end
end

end