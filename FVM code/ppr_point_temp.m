function du = ppr_point_temp(Mesh, PDE, ue)

du = zeros(Mesh.nP, 2);

for P = 1:Mesh.nP
    if Mesh.isbdp(P)
        xc = Mesh.xp(P); yc = Mesh.yp(P);
        du(P,:) = PDE.du(xc, yc);
        continue
    end
    
    nUs = Mesh.P2U{P}; nEs = [];
    for tU = nUs
        nEs = [nEs, Mesh.U2E{tU}];
    end
    nEs = unique(nEs);
        
    tx = Mesh.xe(nEs) - Mesh.xp(P);
    ty = Mesh.ye(nEs) - Mesh.yp(P);

    A = [ones(length(nEs),1), tx, ty, ...
        tx.*tx, tx.*ty, ty.*ty];

%     A = [ones(length(nEs),1), tx, ty, ...
%         tx.*tx, tx.*ty, ty.*ty, ...
%         tx.*tx.*ty, tx.*ty.*ty, tx.*tx.*ty.*ty];
    
    a = A \ ue(nEs);

    du(P,1) = a(2); du(P,2) = a(3);
end
