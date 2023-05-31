function weight = dist_inv_weight(Mesh, ~, ~)
% get interplation weights on mesh points by distance inverse interplation

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
        
        weight{P} = dist_inv / sum(dist_inv);
    end
end

end