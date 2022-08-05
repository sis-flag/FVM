function weight = geo_order2_weight(Mesh, ~, ~)
% get interplation weights on mesh points

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
        
        dist_inv = 1 ./ sqrt((xu-xp).^2 + (yu-yp).^2);
        w0 = dist_inv / sum(dist_inv);
        
        A = zeros(3, length(nU));
        A(1, :) = 1;
        for k = 1:length(nU)
            A(2, k) = xu(k) - xp;
            A(3, k) = yu(k) - yp;
        end
        b = [1; 0; 0];
        
        lambda = (A * A') \ (A * w0 - b);
        w = w0 - A' * lambda;
        
        weight{P} = w / sum(w);
    end
end

end