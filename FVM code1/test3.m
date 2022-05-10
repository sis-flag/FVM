clear
clc

old_path = path;
old_path = path(old_path, './mesh');
old_path = path(old_path, './PDE');
old_path = path(old_path, './interp');

%% Mesh
% Mesh = load_mesh_file('mesh/mesh1_3.txt');

load mesh/mesh11_7
Mesh = arrange_polygonal(coord, U2P);

% Mesh = get_Kershaw_mesh(64, 64, 0.25);

% plot_mesh(Mesh);

%% PDE
PDE = linear_stab1();

u_exact = zeros(Mesh.nU, 1);
for U = 1:Mesh.nU
    xc = Mesh.xc(U); yc = Mesh.yc(U);
    u_exact(U) = PDE.u(xc, yc);
end

%% picard iteration
disp('------picard------')
% u0 = ones(Mesh.nU, 1);
u0 = u_exact;
MaxIt = 1000;

t0 = cputime();
fval = zeros(MaxIt, 1);
for nit = 1:MaxIt
    up = interp_order2_limit(Mesh, PDE, u0);
    [A, F] = mat_FPS(Mesh, PDE, up, u0);
    u1 = A \ F;
    
    if norm(u0-u1, 'inf') < 1e-10
        fval = fval(1:nit);
        break
    elseif any(isnan(u1)) || norm(A * u0 - F, 'inf') > 2e3
        disp('boom!')
        break
    else
        fval(nit) = norm(A * u0 - F, 'inf');
        u0 = u1;
    end
end

disp('Linf error')
disp(max(abs(u_exact - u1)))
disp('Nit')
disp(nit)
disp('time')
disp(cputime() - t0)

%% weighted sparse quasi Newton
disp('------sparse quasi Newton------')
% u0 = ones(Mesh.nU, 1);
u0 = u_exact;
MaxIt = 300;
thres = 0.01;
deepth = 20;

nzInd = cell(Mesh.nU, 1);
for U = 1:Mesh.nU
    for P = Mesh.U2E{U}
        aU = Mesh.E2U{P};
        nzInd{U} = [nzInd{U}, aU];
    end
    nzInd{U} = unique(nzInd{U});
end

t0 = cputime();
fval2 = zeros(MaxIt, 1);
allU = zeros(Mesh.nU, MaxIt);
allP = zeros(Mesh.nU, MaxIt);
for nit = 1:MaxIt
    up = interp_order2_limit(Mesh, PDE, u0);
    [A, F] = mat_FPS(Mesh, PDE, up, u0);
    phi0 = A * u0 - F;
    
    allU(:, nit) = u0;
    allP(:, nit) = phi0;
    
    if nit < deepth+1
        u1 = u0 - A \ phi0;
    else
        dU = diff(allU(:, nit-deepth:nit), 1, 2);
        dP = diff(allP(:, nit-deepth:nit), 1, 2);
        J = approxJ(nzInd, dU, dP);
        
        w = thres / (thres + norm(phi0));
        u1 = u0 - (w*J + (1-w)*A) \ phi0;
    end
    
    if norm(u0-u1, 'inf') < 1e-10
        fval2 = fval2(1:nit);
        break
    elseif any(isnan(u1)) || norm(A * u0 - F, 'inf') > 2e3
        disp('boom!')
        break
    else
        fval2(nit) = norm(phi0, 'inf');
        u0 = u1;
    end
end

disp('Linf error')
disp(max(abs(u_exact - u1)))
disp('Nit')
disp(nit)
disp('time')
disp(cputime() - t0)

%% positive preserving Anderson

%% result
% figure
% plot_func_unit(Mesh, u1)
%
% figure
% u1p = interp_my(Mesh, PDE, u1);
% plot_func_point(Mesh, u1p);
