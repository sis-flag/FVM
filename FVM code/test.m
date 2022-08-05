clear;

old_path = path;
old_path = path(old_path, './mesh');
old_path = path(old_path, './PDE');
old_path = path(old_path, './interp');

%% Mesh
% Mesh = load_mesh_file('mesh/mesh1_2.txt');

load mesh/mesh11_1
Mesh = arrange_polygonal(coord, U2P);

% Mesh = get_Kershaw_mesh(8, 8, 0.1);

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

%% NPS
weight = interp_my_weight(Mesh, PDE);
[A, F] = mat_NPS(Mesh, PDE, weight);
u = A \ F;

disp('NPS Linf error')
disp(max(abs(ueU - u)))

%% FPS
% u0 = zeros(Mesh.nU, 1);
% for nit = 1:300
%     up = interp_order2_limit(Mesh, PDE, u0);
%     [A, F] = mat_FPS(Mesh, PDE, up, u0);
%     u1 = A \ F;
% 
%     if norm(u0 - u1) < 1e-10
%         break
%     else
%         u0 = u1;
%     end
% end
% 
% disp('FPS Linf error')
% disp(max(abs(ueU - u1)))

%% EBS1
[A, F] = mat_EBS1(Mesh, PDE);
u = A \ F;

disp('EBS-1 Linf error')
disp(max(abs(ueE - u)))

%% EBS2
[A, F] = mat_EBS2(Mesh, PDE);
u = A \ F;

disp('EBS-2 Linf error')
disp(max(abs(ueE - u)))

%% plot
% figure
% plot_func_edge(Mesh, u)

% figure
% u1p = interp_average(Mesh, PDE, u1);
% plot_func_point(Mesh, u1p);
