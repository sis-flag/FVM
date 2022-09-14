function du = ppr_continous(Mesh, PDE, ue)

du = zeros(Mesh.nE, 2);

for E = 1:Mesh.nE
    nEs = E; nUs = [];
    while length(nEs) < 6
        for tE = nEs
            nUs = [nUs, Mesh.E2U{tE}];
        end
        nUs = unique(nUs);
        
        for tU = nUs
            nEs = [nEs, Mesh.U2E{tU}];
        end
        nEs = unique(nEs);
    end

%     nUs = Mesh.E2U{E}; nEs = [];
%     if length(nUs) == 1
%         xc = Mesh.xe(E); yc = Mesh.ye(E);
%         du(E,:) = PDE.du(xc, yc);
%         continue
%     end
%     for tU = nUs
%         nEs = [nEs, Mesh.U2E{tU}];
%     end
%     nEs = unique(nEs);
        
    tx = Mesh.xe(nEs) - Mesh.xe(E);
    ty = Mesh.ye(nEs) - Mesh.ye(E);
%     h = (mean(abs(tx)) + mean(abs(ty))) / 2;
%     tx = tx/h; ty = ty/h;

    A = [ones(length(nEs),1), tx, ty, ...
        tx.*tx, tx.*ty, ty.*ty];
%     A = [ones(length(nEs),1), tx, ty, ...
%         tx.*tx, tx.*ty, ty.*ty, ...
%         tx.*tx.*ty, tx.*ty.*ty, tx.*tx.*ty.*ty];
    a = A \ ue(nEs);

    du(E,1) = a(2); du(E,2) = a(3);
end
