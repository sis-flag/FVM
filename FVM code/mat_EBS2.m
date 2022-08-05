function [Amat, F] = mat_EBS2(Mesh, PDE, gamma)

if nargin < 3
    gamma = 1;
end

all_a = cell(Mesh.nU, 1);
for U = 1:Mesh.nU
    xc = Mesh.xc(U); yc = Mesh.yc(U);
    all_a{U} = PDE.a(xc, yc, 0);
end

F = zeros(Mesh.nE, 1);
for E = 1:Mesh.nE
    xe = Mesh.xe(E); ye = Mesh.ye(E);
    
    if length(Mesh.E2U{E}) == 1 % boundary edge
        F(E) = PDE.u(xe, ye);
        
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
for U = 1:Mesh.nU
    nE = Mesh.U2E{U};
    nP = Mesh.U2P{U};
    
    Kinv = 1 / Mesh.area(U);
    xP = Mesh.xp(nP); yP = Mesh.yp(nP);
    xE = Mesh.xe(nE); yE = Mesh.ye(nE);
    xU = Mesh.xc(U); yU = Mesh.yc(U);
    
    n = length(nP);
    R = -eye(n) + diag(ones(n-1,1), 1);
    R(n, 1) = 1;
    
    X = R * [xE, yE];
    N = [xP - xU, yP - yU] * [0,1;-1,0];
    C = eye(n) - X / (X'*X) * X';
    A = Kinv * N * all_a{U} * N' + gamma * (C * C');
    RAR = R' * A * R;
    
    for r = 1:length(nE)
        if length(Mesh.E2U{nE(r)}) > 1
            for c = 1:length(nE)
                nnzA = nnzA + 1;
                iA(nnzA) = nE(r); jA(nnzA) = nE(c);
                vA(nnzA) = RAR(r, c);
            end
        end
    end
end


for E = 1:Mesh.nE
    if length(Mesh.E2U{E}) == 1
        nnzA = nnzA + 1;
        iA(nnzA) = E; jA(nnzA) = E;
        vA(nnzA) = 1;
    end
end

Amat = sparse(iA(1:nnzA), jA(1:nnzA), vA(1:nnzA), ...
    Mesh.nE, Mesh.nE);

end