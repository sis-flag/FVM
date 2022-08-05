function n = norm_edge(Mesh, u, p)

Earea = zeros(Mesh.nE, 1, 'double');
for E = 1:Mesh.nE
    nU = Mesh.E2U{E};
    nP = Mesh.E2P{E};
    
    if length(nU) == 1 % boundary edge
        x1 = Mesh.xc(nU); y1 = Mesh.yc(nU);
        x2 = Mesh.xp(nP(1)); y2 = Mesh.yp(nP(1));
        x3 = Mesh.xp(nP(2)); y3 = Mesh.yp(nP(2));
        
        Earea(E) = abs(0.5 * det([x1-x2, x1-x3; y1-y2, y1-y3]));
        
    else % not boundary edge
        x1 = Mesh.xc(nU(1)); y1 = Mesh.yc(nU(1));
        x3 = Mesh.xc(nU(2)); y3 = Mesh.yc(nU(2));
        x2 = Mesh.xp(nP(1)); y2 = Mesh.yp(nP(1));
        x4 = Mesh.xp(nP(2)); y4 = Mesh.yp(nP(2));
        
        Tarea1 = 0.5 * det([x1-x2, x1-x4; y1-y2, y1-y4]);
        Tarea2 = 0.5 * det([x3-x2, x3-x4; y3-y2, y3-y4]);
        
        Earea(E) = abs(Tarea1) + abs(Tarea2);
    end
end

if p == inf
    n = max(abs(u));
elseif p > 0
    n = sum(Earea .* abs(u).^p)^(1/p);
else
    n = NaN;
end

end