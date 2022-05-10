function up = interp_my(Mesh, PDE, uc)
% get interplation weights on mesh points

all_a = cell(Mesh.nU, 1);
for U = 1:Mesh.nU
    xU = Mesh.xc{U};
    all_a{U} = PDE.a(xU(1), xU(2), uc(U));
end

up = zeros(Mesh.nP, 1);
for P = 1:Mesh.nP
    if Mesh.isbdp(P)
        xP = Mesh.xp{P};
        up(P) = PDE.u(xP(1), xP(2));
    else
        xP = Mesh.xp{P};
        nU = Mesh.P2U{P};
        nE = Mesh.P2E{P};
        
        dist_inv = zeros(length(nU), 1);
        for k = 1:length(nU)
            xU = Mesh.xc{nU(k)};
            dist_inv(k) = 1 / norm(xU - xP);
        end
        w0 = dist_inv / sum(dist_inv);
        
        G = zeros(2*length(nU));
        for k = 1:length(nE)
            nk = Mesh.nv{nE(k)};
            P1 = Mesh.E2P{nE(k)}(1);
            P2 = Mesh.E2P{nE(k)}(2);
            tk = Mesh.xp{P1} - Mesh.xp{P2};
            
            U1 = Mesh.E2U{nE(k)}(1);
            k1 = find(nU == U1);
            U2 = Mesh.E2U{nE(k)}(2);
            k2 = find(nU == U2);
            
            G(2*k-1, 2*k1-1:2*k1) = nk' * all_a{U1};
            G(2*k, 2*k1-1:2*k1) = tk';
            G(2*k-1, 2*k2-1:2*k2) = -nk' * all_a{U2};
            G(2*k, 2*k2-1:2*k2) = -tk';
        end
        [~, ~, V] = svd(G);
        gradU = V(:, end);
        
        A = zeros(3, length(nU));
        A(1, :) = 1;
        for k = 1:length(nU)
            A(2:3, k) = gradU(2*k-1:2*k) .* (Mesh.xc{nU(k)} - Mesh.xp{P});
        end
        b = [1; 0; 0];
        
        lambda = (A * A') \ (A * w0 - b);
        w = w0 - A' * lambda;
        
        if any(w < 0)
            up(P) = sum(uc(nU) .* w0) / sum(w0);
        else
            up(P) = sum(uc(nU) .* w) / sum(w);
        end
    end
end

end