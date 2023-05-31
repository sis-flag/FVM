function flux = flux_ECS2(Mesh, PDE, ue, gamma)

if nargin < 4
    gamma = 1;
end

all_a = cell(Mesh.nU, 1);
for U = 1:Mesh.nU
    xc = Mesh.xc(U); yc = Mesh.yc(U);
    all_a{U} = PDE.a(xc, yc, 0);
end

flux = cell(Mesh.nU, 1);
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

    flux{U} = A * R * ue(nE);
end