function F = ECS1_MF(Mesh, PDE, ue)
% ECS for matrix free solver

flux = flux_ECS1(Mesh, PDE, ue);

F = zeros(Mesh.nE, 1);
for U = 1:Mesh.nU
    nE = Mesh.U2E{U};
    nP = Mesh.U2P{U};
    
    Nn = length(nP);
    for nk = 1:Nn
        E1 = nE(mod(nk-1,Nn)+1);
        E2 = nE(mod(nk,Nn)+1);
        
        F(E1) = F(E1) - flux{U}(nk);
        F(E2) = F(E2) + flux{U}(nk);
    end
end

if PDE.bdtype == 'D'
for E = 1:Mesh.nE
    if length(Mesh.E2U{E}) == 1 % boundary edge
        F(E) = ue(E);
    end
end
end

