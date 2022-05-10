function weight = interp_dist_inv_weight(Mesh, ~)
% get interplation weights on mesh points by distance inverse interplation

weight = cell(Mesh.nP, 1);
for P = 1:Mesh.nP
    if Mesh.isbdp(P)
        weight{P} = [];
    else
        nUs = Mesh.P2U{P};
        dist_inv = zeros(size(nUs));
        for k = 1:length(nUs)
            dist_inv(k) = 1 / norm(Mesh.xc{nUs(k)} - Mesh.xp{P});
        end
        weight{P} = dist_inv / sum(dist_inv);
    end
end

end