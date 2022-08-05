function [Amat, F] = mat_NPS(Mesh, PDE, weight, uc)

if nargin < 4
    uc = zeros(Mesh.nU, 1);
end

F = zeros(Mesh.nU, 1);
for U = 1:Mesh.nU
    xc = Mesh.xc(U); yc = Mesh.yc(U);
    F(U) = Mesh.area(U) * PDE.f(xc, yc, uc(U));
end

all_a = cell(Mesh.nU, 1);
for U = 1:Mesh.nU
    xc = Mesh.xc(U); yc = Mesh.yc(U);
    all_a{U} = PDE.a(xc, yc, uc(U));
end

% Amat = spalloc(Mesh.nU, Mesh.nU, 16*Mesh.nE);
iA = zeros(16*Mesh.nE, 1);
jA = zeros(16*Mesh.nE, 1);
vA = zeros(16*Mesh.nE, 1);
nnzA = 0;
for E = 1:Mesh.nE
    if length(Mesh.E2U{E}) == 1 % boundary edge
        K = Mesh.E2U{E};
        A = Mesh.E2P{E}(1);
        B = Mesh.E2P{E}(2);
        xK = Mesh.xc(K); yK = Mesh.yc(K);
        xA = Mesh.xp(A); yA = Mesh.yp(A);
        xB = Mesh.xp(B); yB = Mesh.yp(B);
        uA = PDE.u(xA, yA);
        uB = PDE.u(xB, yB);
        
        nx = Mesh.nx(E); ny = Mesh.ny(E);
        ane = all_a{K} * [nx; ny];
        xKA = xA - xK; xKB = xB - xK;
        yKA = yA - yK; yKB = yB - yK;
        temp = [xKA, xKB; yKA, yKB] \ ane;
        aKA = temp(1); aKB = temp(2);
        
        % Amat(K, K) += Mesh.len(E) * (aKA + aKB);
        nnzA = nnzA + 1;
        iA(nnzA) = K; jA(nnzA) = K;
        vA(nnzA) = Mesh.len(E) * (aKA + aKB);
        
        F(K) = F(K) + Mesh.len(E) * (aKA*uA + aKB*uB);
        
    elseif length(Mesh.E2U{E}) == 2 % not boundary edge
        K = Mesh.E2U{E}(1);
        L = Mesh.E2U{E}(2);
        A = Mesh.E2P{E}(1);
        B = Mesh.E2P{E}(2);
        xK = Mesh.xc(K); yK = Mesh.yc(K);
        xL = Mesh.xc(L); yL = Mesh.yc(L);
        xA = Mesh.xp(A); yA = Mesh.yp(A);
        xB = Mesh.xp(B); yB = Mesh.yp(B);
        
        nx = Mesh.nx(E); ny = Mesh.ny(E);
        ane = all_a{K} * [nx; ny];
        xKA = xA - xK; xKB = xB - xK;
        yKA = yA - yK; yKB = yB - yK;
        temp = [xKA, xKB; yKA, yKB] \ ane;
        aKA = temp(1); aKB = temp(2);
        
        nx = Mesh.nx(E); ny = Mesh.ny(E);
        ane = all_a{L} * [nx; ny];
        xLA = xA - xL; xLB = xB - xL;
        yLA = yA - yL; yLB = yB - yL;
        temp = [xLA, xLB; yLA, yLB] \ ane;
        aLA = temp(1); aLB = temp(2);
        
        % Amat([K,L], [K,L]) +=  ...
        %     Mesh.len(E)/2 * ...
        %     [aKA+aKB, aLA+aLB; -aKA-aKB , -aLA-aLB];
        
        nnzA = nnzA + 1;
        iA(nnzA) = K; jA(nnzA) = K;
        vA(nnzA) = Mesh.len(E)/2 * (aKA + aKB);
        
        nnzA = nnzA + 1;
        iA(nnzA) = K; jA(nnzA) = L;
        vA(nnzA) = Mesh.len(E)/2 * (aLA + aLB);
        
        nnzA = nnzA + 1;
        iA(nnzA) = L; jA(nnzA) = K;
        vA(nnzA) = - Mesh.len(E)/2 * (aKA + aKB);
        
        nnzA = nnzA + 1;
        iA(nnzA) = L; jA(nnzA) = L;
        vA(nnzA) = - Mesh.len(E)/2 * (aLA + aLB);
        
        if Mesh.isbdp(A)
            uA = PDE.u(xA, yA);
            F(K) = F(K) + Mesh.len(E)/2 * (aKA + aLA) * uA;
            F(L) = F(L) - Mesh.len(E)/2 * (aKA + aLA) * uA;
        else
            Kn = Mesh.P2U{A};
            % Amat([K,L], Kn) = Amat([K,L], Kn) - ...
            %     Mesh.len(E)/2 * (aKA + aLA) * ...
            %     [weight{A}'; -weight{A}'];
            for k = 1:length(Kn)
                nnzA = nnzA + 1;
                iA(nnzA) = K; jA(nnzA) = Kn(k);
                vA(nnzA) = - Mesh.len(E)/2 * (aKA + aLA)...
                    * weight{A}(k);
                
                nnzA = nnzA + 1;
                iA(nnzA) = L; jA(nnzA) = Kn(k);
                vA(nnzA) = Mesh.len(E)/2 * (aKA + aLA) ...
                    * weight{A}(k);
            end
        end
        
        if Mesh.isbdp(B)
            uB = PDE.u(xB, yB);
            F(K) = F(K) + Mesh.len(E)/2 * (aKB + aLB) * uB;
            F(L) = F(L) - Mesh.len(E)/2 * (aKB + aLB) * uB;
        else
            Kn = Mesh.P2U{B};
            % Amat([K,L], Kn) = Amat([K,L], Kn) - ...
            %     Mesh.len(E)/2 * (aKB + aLB) * ...
            %     [weight{B}'; -weight{B}'];
            for k = 1:length(Kn)
                nnzA = nnzA + 1;
                iA(nnzA) = K; jA(nnzA) = Kn(k);
                vA(nnzA) = - Mesh.len(E)/2 * (aKB + aLB)...
                    * weight{B}(k);
                
                nnzA = nnzA + 1;
                iA(nnzA) = L; jA(nnzA) = Kn(k);
                vA(nnzA) = Mesh.len(E)/2 * (aKB + aLB) ...
                    * weight{B}(k);
            end
        end
    end
end

Amat = sparse(iA(1:nnzA), jA(1:nnzA), vA(1:nnzA), ...
    Mesh.nU, Mesh.nU);

end