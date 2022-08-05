function [Amat, F] = mat_FPS(Mesh, PDE, up, uc)

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

    function [a1, a2, P1, P2] = convex_decomp(ane, K)
        xK = Mesh.xc(K); yK = Mesh.yc(K);
        
        nPs = Mesh.U2P{K};
        for k = 1:length(nPs)
            if k == 1
                P1 = nPs(end);
                P2 = nPs(1);
            else
                P1 = nPs(k-1);
                P2 = nPs(k);
            end
            
            x1 = Mesh.xp(P1); y1 = Mesh.yp(P1);
            x2 = Mesh.xp(P2); y2 = Mesh.yp(P2);
            
            xK1 = x1 - xK; xK2 = x2 - xK;
            yK1 = y1 - yK; yK2 = y2 - yK;
            temp = [xK1, xK2; yK1, yK2] \ ane;
            a1 = temp(1); a2 = temp(2);
            
            if a1 > 0 && a2 > 0
                break
            end
        end
    end

% Amat = spalloc(Mesh.nU, Mesh.nU, 16*Mesh.nE);
iA = zeros(16*Mesh.nE, 1);
jA = zeros(16*Mesh.nE, 1);
vA = zeros(16*Mesh.nE, 1);
nnzA = 0;
for E = 1:Mesh.nE
    nx = Mesh.nx(E); ny = Mesh.ny(E);
    
    if length(Mesh.E2U{E}) == 1 % boundary edge
        K = Mesh.E2U{E};
        
        aKne = all_a{K} * [nx; ny];
        [a1, a2, P1, P2] = convex_decomp(aKne, K);
        
        % Amat(K, K) += Mesh.len(E) * (a1 + a2);
        nnzA = nnzA + 1;
        iA(nnzA) = K; jA(nnzA) = K;
        vA(nnzA) = Mesh.len(E) * (a1 + a2);
        
        F(K) = F(K) + Mesh.len(E) * (a1*up(P1) + a2*up(P2));
        
    elseif length(Mesh.E2U{E}) == 2 % not boundary edge
        K = Mesh.E2U{E}(1);
        L = Mesh.E2U{E}(2);
        
        aKne = all_a{K} * [nx; ny];
        aLne = all_a{L} * [nx; ny];
        
        [a1, a2, P1, P2] = convex_decomp(aKne, K);
        [a3, a4, P3, P4] = convex_decomp(-aLne, L);
        
        tK = a1 * up(P1) + a2 * up(P2);
        tL = a3 * up(P3) + a4 * up(P4);
        
        if tK + tL == 0
            muK = 0.5; muL = 0.5;
        else
            muK = tL / (tK + tL);
            muL = tK / (tK + tL);
        end
        
        % Amat([K,L], [K,L]) += ...
        %     Mesh.len{E} * ...
        %     [muK*(a1+a2), -muL*(a3+a4); ...
        %     -muK*(a1+a2), muL*(a3+a4)];
        
        
        nnzA = nnzA + 1;
        iA(nnzA) = K; jA(nnzA) = K;
        vA(nnzA) = Mesh.len(E) * muK * (a1 + a2);
        
        nnzA = nnzA + 1;
        iA(nnzA) = K; jA(nnzA) = L;
        vA(nnzA) = - Mesh.len(E) * muL * (a3 + a4);
        
        nnzA = nnzA + 1;
        iA(nnzA) = L; jA(nnzA) = K;
        vA(nnzA) = - Mesh.len(E) * muK * (a1 + a2);
        
        nnzA = nnzA + 1;
        iA(nnzA) = L; jA(nnzA) = L;
        vA(nnzA) = Mesh.len(E) * muL * (a3 + a4);
    end
end

Amat = sparse(iA(1:nnzA), jA(1:nnzA), vA(1:nnzA), ...
    Mesh.nU, Mesh.nU);

end