function up = interp_limit(Mesh, up, uc)
% monotone limitation. replace exception value

for P = 1:Mesh.nP
    if ~Mesh.isbdp(P)
        
        nU = Mesh.P2U{P};
        
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