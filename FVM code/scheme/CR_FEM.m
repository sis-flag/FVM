function [Amat, F] = CR_FEM(Mesh, PDE)

all_a = cell(Mesh.nU, 1);
for U = 1:Mesh.nU
    xc = Mesh.xc(U); yc = Mesh.yc(U);
    all_a{U} = PDE.a(xc, yc, 0);
end


% Amat = spalloc(Mesh.nU, Mesh.nU, 16*Mesh.nE);
iA = zeros(50*Mesh.nE, 1);
jA = zeros(50*Mesh.nE, 1);
vA = zeros(50*Mesh.nE, 1);
nnzA = 0;
F = zeros(Mesh.nE, 1);
for U = 1:Mesh.nU
    nE = Mesh.U2E{U};
    
    Kinv = 1 / Mesh.area(U);
    xE = Mesh.xe(nE); yE = Mesh.ye(nE);
    xU = Mesh.xc(U); yU = Mesh.yc(U);
    
    temp = 2 * [yE, -xE]' * [0,-1, 1; 1, 0, -1; -1, 1, 0];
    Aloc = Kinv * temp' * all_a{U} * temp;
    
    Floc = Mesh.area(U)/3*[1;1;1] * PDE.f(xU, yU);
    
    for r = 1:length(nE)
        if length(Mesh.E2U{nE(r)}) > 1
            for c = 1:length(nE)
                nnzA = nnzA + 1;
                iA(nnzA) = nE(r); jA(nnzA) = nE(c);
                vA(nnzA) = Aloc(r, c);
            end
        end
    end
    
    F(nE) = F(nE) + Floc;
end


for E = 1:Mesh.nE
    % boundary edge
    if length(Mesh.E2U{E}) == 1
        nnzA = nnzA + 1;
        iA(nnzA) = E; jA(nnzA) = E;
        vA(nnzA) = 1;
        
        xe = Mesh.xe(E); ye = Mesh.ye(E);
        F(E) = PDE.u(xe, ye);
    end
end

Amat = sparse(iA(1:nnzA), jA(1:nnzA), vA(1:nnzA), ...
    Mesh.nE, Mesh.nE);


end