function flux = flux_ECS1(Mesh, PDE, ue)

all_a = cell(Mesh.nU, 1);
for U = 1:Mesh.nU
    xc = Mesh.xc(U); yc = Mesh.yc(U);
    all_a{U} = PDE.a(xc, yc, 0);
end

flux1 = cell(Mesh.nU, 1);
flux2 = cell(Mesh.nU, 1);
for U = 1:Mesh.nU
    nE = Mesh.U2E{U};
    nP = Mesh.U2P{U};
    
    Nn = length(nP);
    for nk = 1:Nn
        P = nP(nk);
        
        E0 = nE(mod(nk-2,Nn)+1);
        E1 = nE(mod(nk-1,Nn)+1);
        E2 = nE(mod(nk,Nn)+1);
        E3 = nE(mod(nk+1,Nn)+1);
        
        xU = Mesh.xc(U); yU = Mesh.yc(U);
        xP = Mesh.xp(P); yP = Mesh.yp(P);
        x0 = Mesh.xe(E0); y0 = Mesh.ye(E0);
        x1 = Mesh.xe(E1); y1 = Mesh.ye(E1);
        x2 = Mesh.xe(E2); y2 = Mesh.ye(E2);
        x3 = Mesh.xe(E3); y3 = Mesh.ye(E3);
        
        ane = all_a{U} * [yU-yP; xP-xU];
        
        temp = [x0-x1, x2-x1; y0-y1, y2-y1] \ ane;
        a10 = temp(1); a12 = temp(2);
        
        temp = [x1-x2, x3-x2; y1-y2, y3-y2] \ ane;
        a21 = temp(1); a23 = temp(2);

        flux1{U}(nk) = a10 * (ue(E1) - ue(E0)) + a12 * (ue(E1) - ue(E2));
        flux2{U}(nk) = a21 * (ue(E2) - ue(E1)) + a23 * (ue(E2) - ue(E3));
    end
end

flux = cell(Mesh.nU, 1);
for U = 1:Mesh.nU
    flux{U} = -(flux1{U} + flux2{U}) / 2;
end