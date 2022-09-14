function du = ppr_continous_point(Mesh, PDE, ue)

du = zeros(Mesh.nP, 2);

for P = 1:Mesh.nP
    xc = Mesh.xp(P); yc = Mesh.yp(P);
    if xc < 1/4 || xc > 3/4 || yc < 1/4 || yc > 3/4
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
%     h = (mean(abs(tx)) + mean(abs(ty))) / 2;
%     tx = tx/h; ty = ty/h;

    A = [ones(length(nEs),1), tx, ty, ...
        tx.*tx, tx.*ty, ty.*ty];
%     A = [ones(length(nEs),1), tx, ty, ...
%         tx.*tx, tx.*ty, ty.*ty, ...
%         tx.*tx.*ty, tx.*ty.*ty, tx.*tx.*ty.*ty];
    a = A \ ue(nEs);

    du(P,1) = a(2); du(P,2) = a(3);
end
