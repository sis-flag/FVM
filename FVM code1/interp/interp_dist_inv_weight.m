function weight = interp_dist_inv_weight(Mesh, ~, ~)
% get interplation weights on mesh points by distance inverse interplation

weight = cell(Mesh.nP, 1);
for P = 1:Mesh.nP
    if Mesh.isbdp(P)
        weight{P} = [];
    else
        nU = Mesh.P2U{P};
        xp = Mesh.xp(P); yp = Mesh.yp(P);
        
        dist_inv = zeros(length(nU), 1);
        for k = 1:length(nU)
            xu = Mesh.xc(nU(k)); yu = Mesh.yc(nU(k));
            dist_inv(k) = 1 / norm([xu-xp; yu-yp]);
        end
        
        weight{P} = dist_inv / sum(dist_inv);
    end
end

end