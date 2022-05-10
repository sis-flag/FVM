function weight = interp_my_weight(Mesh, PDE)
% get interplation weights on mesh points

% TODO: discontinous K

weight = cell(Mesh.nP, 1);
for P = 1:Mesh.nP
    if Mesh.isbdp(P)
        weight{P} = [];
    else
        nUs = Mesh.P2U{P};
        
        A = zeros(3, length(nUs));
        A(1, :) = 1;
        for k = 1:length(nUs)
            A(2:3, k) = Mesh.xc{nUs(k)} - Mesh.xp{P};
        end
        b = [1; 0; 0];
        w0 = ones(1, length(nUs)) / length(nUs);
        
        lambda = (A * A') \ (A * w0 - b);
        w = w_0 - A' * lambda;
        
        if any(w < 0)
            dist_inv = zeros(size(nUs));
            for k = 1:length(nUs)
                dist_inv(k) = 1 / norm(Mesh.xc{nUs(k)} - Mesh.xp{P});
            end
            weight{P} = dist_inv / sum(dist_inv);
        else
            weight{P} = w';
        end
    end
end

end

end