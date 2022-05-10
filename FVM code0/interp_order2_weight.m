function weight = interp_order2_weight(Mesh, PDE)
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
        w0 = ones(1, length(nUs))' / length(nUs);
        
        lambda = (A * A') \ (A * w0 - b);
        w = w0 - A' * lambda;
        
        if any(w < 0)
            warning('not monotone')
        end
        
        weight{P} = w';
    end
end

end