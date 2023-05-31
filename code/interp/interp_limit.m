function up = interp_limit(Mesh, up, uc)
% monotone limitation. replace exception value

for P = 1:Mesh.nP
    if ~Mesh.isbdp(P)
        nU = Mesh.P2U{P};
        
        % monotone limit
        maxu = max(uc(nU));
        minu = min(uc(nU));
        if up(P) > maxu
            up(P) = maxu;
        elseif up(P) < minu
            up(P) = minu;
        end
    end
end

end