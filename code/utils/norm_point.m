function n = norm_point(Mesh, u, p)

Parea = zeros(Mesh.nP, 1);
for P = 1:Mesh.nP
    xp = Mesh.xp(P); yp = Mesh.yp(P);
    
    nEs = Mesh.P2E{P};
    nUs = Mesh.P2U{P};

    for U = nUs
        tnE = intersect(nEs, Mesh.U2E{U});

        xe1 = Mesh.xe(tnE(1)); ye1 = Mesh.ye(tnE(1));
        xe2 = Mesh.xe(tnE(2)); ye2 = Mesh.ye(tnE(2));
        xc = Mesh.xc(U); yc = Mesh.yc(U);

        Tarea1 = det([xe1-xp, xe1-xc; ye1-yp, ye1-yc]);
        Tarea2 = det([xe2-xp, xe2-xc; ye2-yp, ye2-yc]);

        Parea(P) = Parea(P) + 0.5 * abs(Tarea1) + 0.5 * abs(Tarea2);
    end
end

if p == inf
    n = max(abs(u));
elseif p > 0
    n = sum(Parea .* abs(u).^p)^(1/p);
else
    error('p is not valid')
end

end