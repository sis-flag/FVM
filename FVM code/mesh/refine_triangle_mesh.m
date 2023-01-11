function Mesh = refine_triangle_mesh(Mesh0)

N = Mesh0.nP;
coord = [Mesh0.xp, Mesh0.yp; Mesh0.xe, Mesh0.ye]';

U2P = cell(4*Mesh0.nU, 1);
for U = 1:Mesh0.nU
    nE = Mesh0.U2E{U} + N;
    nP = Mesh0.U2P{U};
    
    if length(nE) ~= 3 || length(nP) ~= 3
        error('not triangle!')
    end
    
    U2P{4*U-3} = [nP(1), nE(2), nE(1)];
    U2P{4*U-2} = [nE(2), nP(2), nE(3)];
    U2P{4*U-1} = [nE(1), nE(3), nP(3)];
    U2P{4*U} = nE; 
end

Mesh = arrange_polygonal(coord, U2P);

end