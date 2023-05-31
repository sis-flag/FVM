clear
addpath("mesh", "PDE", "scheme", "interp", "utils")

% Mesh = get_sin_tri_mesh(4, 4, 0);
% PDE = linear_stab10();
% U_test = 6;

delta = 1;
% Mesh = get_sin_tri_mesh(8, 8, 0);
Mesh = load_mesh_file('mesh/mesh1_2.txt');
PDE = struct();
PDE.u = @(x,y) sin(pi*x) * sin(pi*y);
PDE.a = @(~,~,~) [delta, 0; 0, 1];
PDE.du = @(x,y) pi * [cos(pi*x)*sin(pi*y); cos(pi*y)*sin(pi*x)];
PDE.f = @(x,y,~) sin(pi*x) * sin(pi*y) * (delta+1)*pi*pi;
PDE.bdtype = "D";
U_test = 2;

% delta = 1000;
% Mesh = get_sin_tri_mesh(2, 2, 0);
% PDE = struct();
% PDE.u = @(x,y) x*(1-x)*y*(1-y);
% PDE.a = @(~,~,~) [delta, 0; 0, 1];
% PDE.du = @(x,y)  [y*(2*x-1)*(y-1); x*(2*y-1)*(x-1)];
% PDE.f = @(x,y,~) -2*x*(x-1) - 2*delta*y*(y-1);
% PDE.bdtype = "D";

% plot_mesh(Mesh)

u_exact = zeros(Mesh.nE, 1);
for E = 1:Mesh.nE
    xc = Mesh.xe(E); yc = Mesh.ye(E);
    u_exact(E) = PDE.u(xc, yc);
end

du_exact = cell(Mesh.nU, 1);
t = sym("t");
for U = U_test
    nP = Mesh.U2P{U};
    
    Nn = length(nP);
    for nk = 1:Nn
        P = nP(nk);
        
        xU = Mesh.xc(U); yU = Mesh.yc(U);
        xP = Mesh.xp(P); yP = Mesh.yp(P);

        du_sym = PDE.du(t*xP+(1-t)*xU, t*yP+(1-t)*yU);
        D_sym = int(du_sym, t, 0, 1);
        
        du_exact{U}(:, nk) = double(D_sym);
    end
end

% flux_exact = cell(Mesh.nU, 1);
% for U = 1:U_test
%     nP = Mesh.U2P{U};
%     
%     Nn = length(nP);
%     for nk = 1:Nn
%         P = nP(nk);
%         
%         xU = Mesh.xc(U); yU = Mesh.yc(U);
%         xP = Mesh.xp(P); yP = Mesh.yp(P);
% 
%         flux_exact{U}(nk) = [yU-yP, xP-xU] * ...
%             PDE.a(0,0) * du_exact{U}(:, nk);
%     end
% end

% flux_exact = cell(Mesh.nU, 1);
% t = sym("t");
% for U = 1:Mesh.nU
%     nP = Mesh.U2P{U};
%     
%     Nn = length(nP);
%     for nk = 1:Nn
%         P = nP(nk);
%         
%         xU = Mesh.xc(U); yU = Mesh.yc(U);
%         xP = Mesh.xp(P); yP = Mesh.yp(P);
% 
%         du_sym = PDE.du(t*xP+(1-t)*xU, t*yP+(1-t)*yU);
%         a_sym = PDE.a(t*xP+(1-t)*xU, t*yP+(1-t)*yU);
%         temp = [yU-yP, xP-xU] * a_sym * du_sym;
%         F_sym = int(temp, t, 0, 1);
%         
%         flux_exact{U}(nk) = double(F_sym);
%     end
% end
% flux_exact = cell2mat(flux_exact);

% F_exact = zeros(Mesh.nE, 1);
% syms xi eta
% for E = 1:Mesh.nE
%     xe = Mesh.xe(E); ye = Mesh.ye(E);
%     
%     if length(Mesh.E2U{E}) == 1 % boundary edge
%         if PDE.bdtype == 'D' % Dirichlet boundary
%             F_exact(E) = PDE.u(xe, ye);
%         end
%     else % not boundary edge (must parallegram)
%         nU = Mesh.E2U{E};
%         x1 = Mesh.xc(nU(1)); y1 = Mesh.yc(nU(1));
%         x3 = Mesh.xc(nU(2)); y3 = Mesh.yc(nU(2));
%         
%         nE = Mesh.E2P{E};
%         x2 = Mesh.xp(nE(1)); y2 = Mesh.yp(nE(1));
%         x4 = Mesh.xp(nE(2)); y4 = Mesh.yp(nE(2));
% 
%         tx = [x2-x1, x4-x1; y2-y1, y4-y1] * [xi; eta] + [x1; y1];
%         J = det([x2-x1, x4-x1; y2-y1, y4-y1]);
% 
%         temp = PDE.f(tx(1), tx(2), 0);
%         temp = int(temp, xi, 0, 1);
%         temp = int(temp, eta, 0, 1);
% 
%         F_exact(E) = J * double(temp);
%     end
% end

% flux_test = flux_ECS1(Mesh, PDE, u_exact);
% flux_test = cell2mat(flux_test);

du_test = cell(Mesh.nU, 1);
for U = U_test
    nE = Mesh.U2E{U}; nP = Mesh.U2P{U}; uu = u_exact(nE);
    x1 = Mesh.xe(nE(1)); y1 = Mesh.ye(nE(1));
    x2 = Mesh.xe(nE(2)); y2 = Mesh.ye(nE(2));
    x3 = Mesh.xe(nE(3)); y3 = Mesh.ye(nE(3));
    dduu = [x1-x2, y1-y2; x1-x3, y1-y3] \ [uu(1) - uu(2); uu(1) - uu(3)];

    du_test{U} = dduu;
end

[A, F] = ECS1(Mesh, PDE);
u_ECS = A \ F;


flux_ECS = flux_ECS1(Mesh, PDE, u_ECS);
flux_ECS = cell2mat(flux_ECS);

% for U = 1:Mesh.nU
% disp([U, flux_test{U} - flux_ECS{U}])
% end

% figure
% pcolor_func_edge(Mesh, A * u_exact - F)
% 
% figure
% pcolor_func_edge(Mesh, u_ECS - u_exact)


% % WTF ???
% flux_test{U_test}
% flux_exact{U_test}
% 
% nP = Mesh.U2P{U_test}
% nE = Mesh.U2E{U_test}
% 
% x = Mesh.xc(U_test);
% y = Mesh.yc(U_test);
% a = PDE.a(x,y)
% 
% u_exact(nE)