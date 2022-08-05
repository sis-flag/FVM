function weight = order2_weight(Mesh, PDE, uc)
% get interplation weights on mesh points

if nargin < 3
    uc = zeros(Mesh.nU, 1);
end

all_a = cell(Mesh.nU, 1);
for U = 1:Mesh.nU
    xu = Mesh.xc(U);
    yu = Mesh.yc(U);
    all_a{U} = PDE.a(xu, yu, uc(U));
end

weight = cell(Mesh.nP, 1);
for P = 1:Mesh.nP
    if Mesh.isbdp(P)
        weight{P} = [];
    else
        xp = Mesh.xp(P);
        yp = Mesh.yp(P);
        
        nU = Mesh.P2U{P};
        xu = Mesh.xc(nU);
        yu = Mesh.yc(nU);
        
        nE = Mesh.P2E{P};
        
        dist_inv = 1 ./ sqrt((xu-xp).^2 + (yu-yp).^2);
        w0 = dist_inv / sum(dist_inv);
        
        G = zeros(2*length(nU));
        for k = 1:length(nE)
            nk = [Mesh.nx(nE(k)); Mesh.ny(nE(k))];
            P1 = Mesh.E2P{nE(k)}(1);
            P2 = Mesh.E2P{nE(k)}(2);
            xp1 = Mesh.xp(P1);
            xp2 = Mesh.xp(P2);
            yp1 = Mesh.yp(P1);
            yp2 = Mesh.yp(P2);
            tk = [xp1-xp2; yp1-yp2];
            
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
            A(2, k) = gradU(2*k-1) * (xu(k) - xp);
            A(3, k) = gradU(2*k) * (yu(k) - yp);
        end
        b = [1; 0; 0];
        
        lambda = (A * A') \ (A * w0 - b);
        w = w0 - A' * lambda;
        
        weight{P} = w / sum(w);
    end
end

end