function uc = solve_NPS(Mesh, PDE, weight)

if PDE.bdtype == 'N'
    F = zeros(Mesh.nU+1, 1);
else
    F = zeros(Mesh.nU, 1);
end
for U = 1:Mesh.nU
    xcu = Mesh.xc{U};
    F(U) = Mesh.area{U} * PDE.f(xcu(1), xcu(2));
end

all_a = cell(Mesh.nU, 1);
for U = 1:Mesh.nU
    xcu = Mesh.xc{U};
    all_a{U} = PDE.a(xcu(1), xcu(2));
end

if PDE.bdtype == 'N'
    M = spalloc(Mesh.nU+1, Mesh.nU, 16*Mesh.nE);
else
    M = spalloc(Mesh.nU, Mesh.nU, 16*Mesh.nE);
end
for E = 1:Mesh.nE
    if length(Mesh.E2U{E}) == 1 % boundary edge
        if PDE.bdtype == 'D' % Dirichlet boundary
            K = Mesh.E2U{E};
            A = Mesh.E2P{E}(1);
            B = Mesh.E2P{E}(2);
            xK = Mesh.xc{K};
            xA = Mesh.xp{A};
            xB = Mesh.xp{B};
            uA = PDE.u(xA(1), xA(2));
            uB = PDE.u(xB(1), xB(2));
            
            KA = xA - xK; KB = xB - xK;
            ane = all_a{K} * Mesh.nv{E};
            temp = [KA, KB] \ ane;
            aKA = temp(1); aKB = temp(2);
            
            M(K, K) = M(K, K) + Mesh.len{E} * (aKA + aKB);
            F(K) = F(K) + Mesh.len{E} * (aKA*uA + aKB*uB);
        elseif PDE.bdtype == 'N' % Neumann boundary
            xe = Mesh.xe{E};
            F(K) = F(K) + Mesh.len{E} * PDE.g(xe(1), xe(2));
        end
    elseif length(Mesh.E2U{E}) == 2 % not boundary edge
        K = Mesh.E2U{E}(1);
        L = Mesh.E2U{E}(2);
        A = Mesh.E2P{E}(1);
        B = Mesh.E2P{E}(2);
        xK = Mesh.xc{K};
        xL = Mesh.xc{L};
        xA = Mesh.xp{A};
        xB = Mesh.xp{B};
        
        KA = xA - xK; KB = xB - xK;
        ane = all_a{K} * Mesh.nv{E};
        temp = [KA, KB] \ ane;
        aKA = temp(1); aKB = temp(2);
        
        LA = xA - xL; LB = xB - xL;
        ane = all_a{L} * Mesh.nv{E};
        temp = [LA, LB] \ ane;
        aLA = temp(1); aLB = temp(2);
        
        M([K,L], [K,L]) = M([K,L], [K,L]) + ...
            Mesh.len{E}/2 * ...
            [aKA+aKB, aLA+aLB; -aKA-aKB , -aLA-aLB];
        
        if Mesh.isbdp(A)
            uA = PDE.u(xA(1), xA(2));
            F([K,L]) = F([K,L]) + ...
                Mesh.len{E}/2 * (aKA + aLA) * ...
                [uA; -uA];
        else
            Kn = Mesh.P2U{A};
            M([K,L], Kn) = M([K,L], Kn) - ...
                Mesh.len{E}/2 * (aKA + aLA) * ...
                [weight{A}; -weight{A}];
        end
        
        if Mesh.isbdp(B)
            uB = PDE.u(xB(1), xB(2));
            F([K,L]) = F([K,L]) + ...
                Mesh.len{E}/2 * (aKB + aLB) * ...
                [uB; -uB];
        else
            Kn = Mesh.P2U{B};
            M([K,L], Kn) = M([K,L], Kn) - ...
                Mesh.len{E}/2 * (aKB + aLB) * ...
                [weight{B}; -weight{B}];
        end
    end
end

if PDE.bdtype == 'N'
    for U = 1:Mesh.nU
        M(end, U) = Mesh.area{U};
    end
end
    
uc = M \ F;
    
end