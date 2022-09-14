clear;

old_path = path;
old_path = path(old_path, './mesh');
old_path = path(old_path, './PDE');
old_path = path(old_path, './interp');
old_path = path(old_path, './scheme');
old_path = path(old_path, './utils');

%% Mesh
% Mesh = load_mesh_file('mesh/mesh1_1.txt');

% load mesh/mesh11_3
% Mesh = arrange_polygonal(coord, U2P);

Mesh = get_Kershaw_mesh(16, 16, 0.4);

plot_mesh(Mesh);

%% PDE
PDE = linear_stab1();

ueU = zeros(Mesh.nU, 1);
for U = 1:Mesh.nU
    xc = Mesh.xc(U); yc = Mesh.yc(U);
    ueU(U) = PDE.u(xc, yc);
end

ueE = zeros(Mesh.nU, 1);
for E = 1:Mesh.nE
    xc = Mesh.xe(E); yc = Mesh.ye(E);
    ueE(E) = PDE.u(xc, yc);
end

dueE = zeros(Mesh.nU, 2);
for E = 1:Mesh.nE
    xc = Mesh.xe(E); yc = Mesh.ye(E);
    dueE(E, :) = PDE.du(xc, yc);
end

%% NPS
% weight = geo_order2_weight(Mesh, PDE);
% [A, F] = NPS(Mesh, PDE, weight);
% u = A \ F;
% 
% disp('NPS Linf error')
% disp(max(abs(ueU - u)))

%% FPS
% u0 = zeros(Mesh.nU, 1);
% for nit = 1:300
%     weight = order2_weight(Mesh, PDE, u0);
%     up = interp_by_weight(Mesh, PDE, u0, weight);
%     up = interp_limit(Mesh, up, u0);
%     [A, F] = FPS(Mesh, PDE, up, u0);
%     u1 = A \ F;
% 
%     if norm(u0-u1) < 1e-10
%         break
%     else
%         u0 = u1;
%     end
% end
% 
% disp('FPS Linf error')
% disp(max(abs(ueU - u1)))
% disp('FPS iteration')
% disp(nit)

%% EBS1
% [A, F] = ECS1(Mesh, PDE);
% u = A \ F;
% 
% disp('ECS-I Linf error')
% disp(max(abs(ueE - u)))

%% EBS2
% [A, F] = ECS2(Mesh, PDE);
% u = A \ F;
% 
% disp('ECS-II Linf error')
% disp(max(abs(ueE - u)))

%% PPR
du = ppr_continous(Mesh, PDE, ueE);

disp('PPR error')
disp(max(max(abs(dueE - du))))

%% plot
% figure
% plot_func_edge(Mesh, u)

% figure
% u1p = interp_average(Mesh, PDE, u1);
% plot_func_point(Mesh, u1p);
