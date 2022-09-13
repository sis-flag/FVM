function n = norm_point(Mesh, u, p)

Parea = zeros(Mesh.nP, 1, 'double');
for P = 1:Mesh.nP
    nU = Mesh.P2U{P};
    nE = Mesh.P2E{P};
    
    sx = zeros(2*length(nE),1);
    sy = zeros(2*length(nE),1);
    if Mesh.isbdp(P) % boundary point
        % TODO
    else % not boundary point
        % TODO
    end
    
    Parea(P) = get_area(sx, sy);
end

if p == inf
    n = max(abs(u));
elseif p > 0
    n = sum(Parea .* abs(u).^p)^(1/p);
else
    error('p is not valid')
end

end