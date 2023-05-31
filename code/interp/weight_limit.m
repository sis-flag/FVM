function weight = weight_limit(Mesh, weight)
% monotone correction. replace negtive weight by dist_inv weight

for P = 1:Mesh.nP
    if Mesh.isbdp(P)
        weight{P} = [];
    elseif any(weight{P} < 0)
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