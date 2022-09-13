function n = norm_edge(Mesh, u, p)

Earea = zeros(Mesh.nE, 1, 'double');
for E = 1:Mesh.nE
    nU = Mesh.E2U{E};
    nP = Mesh.E2P{E};
    
    if length(nU) == 1 % boundary edge
        sx = zeros(3,1); sy = zeros(3,1);
        sx(1) = Mesh.xc(nU); sy(1) = Mesh.yc(nU);
        sx(2:3) = Mesh.xp(nP); sy(2:3) = Mesh.yp(nP);

    else % not boundary edge
        sx = zeros(4,1); sy = zeros(4,1);
        sx([1,3]) = Mesh.xc(nU); sy([1,3]) = Mesh.yc(nU);
        sx([2,4]) = Mesh.xp(nP); sy([2,4]) = Mesh.yp(nP);
    end
    
    Earea(E) = get_area(sx, sy);
end

if p == inf
    n = max(abs(u));
elseif p > 0
    n = sum(Earea .* abs(u).^p)^(1/p);
else
    error('p is not valid')
end

end