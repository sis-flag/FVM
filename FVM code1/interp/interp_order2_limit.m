function up = interp_order2_limit(Mesh, PDE, uc)
% get interplation weights on mesh points

up = zeros(Mesh.nP, 1);
for P = 1:Mesh.nP
    if Mesh.isbdp(P)
        xp = Mesh.xp(P); yp = Mesh.yp(P);
        up(P) = PDE.u(xp, yp);
    else
        nU = Mesh.P2U{P};
        xp = Mesh.xp(P); yp = Mesh.yp(P);
        
        w0 = ones(length(nU), 1) / length(nU);
        
        A = zeros(3, length(nU));
        A(1, :) = 1;
        for k = 1:length(nU)
            xu = Mesh.xc(nU(k));
            yu = Mesh.yc(nU(k));
            A(2, k) = xu - xp;
            A(3, k) = yu - yp;
        end
        
        b = [1; 0; 0];
        
        lambda = (A * A') \ (A * w0 - b);
        w = w0 - A' * lambda;
        
        up(P) = sum(uc(nU) .* w) / sum(w);
        
        % monotone limit
        Mu = max(uc(nU));
        mu = min(uc(nU));
        if up(P) > Mu
            up(P) = Mu;
        elseif up(P) < mu
            up(P) = mu;
        end
    end
end

end