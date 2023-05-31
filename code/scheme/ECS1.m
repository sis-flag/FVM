function [Amat, F] = ECS1(Mesh, PDE)

all_a = cell(Mesh.nU, 1);
for U = 1:Mesh.nU
    xc = Mesh.xc(U); yc = Mesh.yc(U);
    all_a{U} = PDE.a(xc, yc, 0);
end

F = zeros(Mesh.nE, 1);
for E = 1:Mesh.nE
    xe = Mesh.xe(E); ye = Mesh.ye(E);
    
    if length(Mesh.E2U{E}) == 1 % boundary edge
        if PDE.bdtype == 'D' % Dirichlet boundary
            F(E) = PDE.u(xe, ye);

        elseif PDE.bdtype == 'N' % Neumann boundary
            nU = Mesh.E2U{E};
            x1 = Mesh.xc(nU); y1 = Mesh.yc(nU);

            nP = Mesh.E2P{E};
            xe = Mesh.xe(E); ye = Mesh.ye(E);
            ne = [Mesh.nx(E); Mesh.ny(E)];
            x2 = Mesh.xp(nP(1)); y2 = Mesh.yp(nP(1));
            x4 = Mesh.xp(nP(2)); y4 = Mesh.yp(nP(2));

            xt1 = (x1 + x2 + x4) / 3;
            yt1 = (y1 + y2 + y4) / 3;
            Tarea1 = det([x1-x2, x1-x4; y1-y2, y1-y4]);
            F1 = 0.5 * abs(Tarea1) * PDE.f(xt1, yt1, 0);

            F2 = - Mesh.len(E) * (PDE.du(xe, ye)' * PDE.a(xe, ye) * ne);

            F(E) = F1 + F2;
        end
    else % not boundary edge
        nU = Mesh.E2U{E};
        x1 = Mesh.xc(nU(1)); y1 = Mesh.yc(nU(1));
        x3 = Mesh.xc(nU(2)); y3 = Mesh.yc(nU(2));
        
        nP = Mesh.E2P{E};
        x2 = Mesh.xp(nP(1)); y2 = Mesh.yp(nP(1));
        x4 = Mesh.xp(nP(2)); y4 = Mesh.yp(nP(2));
        
        xt1 = (x1 + x2 + x4) / 3;
        yt1 = (y1 + y2 + y4) / 3;
        Tarea1 = det([x1-x2, x1-x4; y1-y2, y1-y4]);
        F1 = 0.5 * abs(Tarea1) * PDE.f(xt1, yt1, 0);
        
        xt2 = (x3 + x2 + x4) / 3;
        yt2 = (y3 + y2 + y4) / 3;
        Tarea2 = det([x3-x2, x3-x4; y3-y2, y3-y4]);
        F2 = 0.5 * abs(Tarea2) * PDE.f(xt2, yt2, 0);
        
        F(E) = F1 + F2;
    end
end

% Amat = spalloc(Mesh.nU, Mesh.nU, 16*Mesh.nE);
iA = zeros(50*Mesh.nE, 1);
jA = zeros(50*Mesh.nE, 1);
vA = zeros(50*Mesh.nE, 1);
nnzA = 0;

    function Amat_add(i, j, v)
        nnzA = nnzA + 1;
        iA(nnzA) = i;
        jA(nnzA) = j;
        vA(nnzA) = v;
    end

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
        
        % if not boundary
        % Amat([E1, E2], [E0, E1, E2]) += ...
        %     0.5 * [-a10, a10+a12 ,-a12; ...
        %     a10, -a10-a12 , a12];
        % Amat([E1, E2], [E1, E2, E3]) += ...
        %     0.5 * [-a21, a21+a23 ,-a23; ...
        %     a21, -a21-a23 , a23];
        
        if length(Mesh.E2U{E1}) == 2 || PDE.bdtype == 'N'
            Amat_add(E1, E0, -0.5 * a10);
            Amat_add(E1, E1, 0.5 * (a10+a12));
            Amat_add(E1, E2, -0.5 * a12);
            
            Amat_add(E1, E1, -0.5 * a21);
            Amat_add(E1, E2, 0.5 * (a21+a23));
            Amat_add(E1, E3, -0.5 * a23);
        end
        
        if length(Mesh.E2U{E2}) == 2 || PDE.bdtype == 'N'
            Amat_add(E2, E0, 0.5 * a10);
            Amat_add(E2, E1, -0.5 * (a10+a12));
            Amat_add(E2, E2, 0.5 * a12);
            
            Amat_add(E2, E1, 0.5 * a21);
            Amat_add(E2, E2, -0.5 * (a21+a23));
            Amat_add(E2, E3, 0.5 * a23);
        end
    end
end

if PDE.bdtype == 'D'
for E = 1:Mesh.nE
    if length(Mesh.E2U{E}) == 1
        Amat_add(E, E, 1);
    end
end
end

Amat = sparse(iA(1:nnzA), jA(1:nnzA), vA(1:nnzA), ...
    Mesh.nE, Mesh.nE);

end