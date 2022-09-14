clear;

old_path = path;
old_path = path(old_path, './mesh');
old_path = path(old_path, './PDE');
old_path = path(old_path, './interp');
old_path = path(old_path, './scheme');
old_path = path(old_path, './utils');

PDE = linear_stab1();

% p = 0.0;
% all_Mesh = cell(5,1);
% all_Mesh{1} = get_rand_mesh(8, 8, p);
% all_Mesh{2} = get_rand_mesh(16, 16, p);
% all_Mesh{3} = get_rand_mesh(32, 32, p);
% all_Mesh{4} = get_rand_mesh(64, 64, p);
% all_Mesh{5} = get_rand_mesh(128, 128, p);

% all_Mesh = cell(7,1);
% load mesh/mesh11_1
% all_Mesh{1} = arrange_polygonal(coord, U2P);
% load mesh/mesh11_2
% all_Mesh{2} = arrange_polygonal(coord, U2P);
% load mesh/mesh11_3
% all_Mesh{3} = arrange_polygonal(coord, U2P);
% load mesh/mesh11_4
% all_Mesh{4} = arrange_polygonal(coord, U2P);
% load mesh/mesh11_5
% all_Mesh{5} = arrange_polygonal(coord, U2P);
% load mesh/mesh11_6
% all_Mesh{6} = arrange_polygonal(coord, U2P);
% load mesh/mesh11_7
% all_Mesh{7} = arrange_polygonal(coord, U2P);

% all_Mesh = cell(6,1);
% all_Mesh{1} = load_mesh_file('mesh/mesh4_1_1.txt');
% all_Mesh{2} = load_mesh_file('mesh/mesh4_1_2.txt');
% all_Mesh{3} = load_mesh_file('mesh/mesh4_1_3.txt');
% all_Mesh{4} = load_mesh_file('mesh/mesh4_1_4.txt');
% all_Mesh{5} = load_mesh_file('mesh/mesh4_1_5.txt');
% all_Mesh{6} = load_mesh_file('mesh/mesh4_1_6.txt');

% all_Mesh = cell(5,1);
% all_Mesh{1} = load_mesh_file('mesh/mesh1_1.txt');
% all_Mesh{2} = load_mesh_file('mesh/mesh1_2.txt');
% all_Mesh{3} = load_mesh_file('mesh/mesh1_3.txt');
% all_Mesh{4} = load_mesh_file('mesh/mesh1_4.txt');
% all_Mesh{5} = load_mesh_file('mesh/mesh1_5.txt');

p = 0.4;
all_Mesh = cell(6,1);
all_Mesh{1} = get_Kershaw_mesh(4, 4, p);
all_Mesh{2} = get_Kershaw_mesh(8, 8, p);
all_Mesh{3} = get_Kershaw_mesh(16, 16, p);
all_Mesh{4} = get_Kershaw_mesh(32, 32, p);
all_Mesh{5} = get_Kershaw_mesh(64, 64, p);
all_Mesh{6} = get_Kershaw_mesh(128, 128, p);

% p = 0.4;
% all_Mesh = cell(6,1);
% all_Mesh{1} = get_sin_mesh(4, 4, p);
% all_Mesh{2} = get_sin_mesh(8, 8, p);
% all_Mesh{3} = get_sin_mesh(16, 16, p);
% all_Mesh{4} = get_sin_mesh(32, 32, p);
% all_Mesh{5} = get_sin_mesh(64, 64, p);
% all_Mesh{6} = get_sin_mesh(128, 128, p);

%% NPS
% fprintf('\nNPS:\n')
% fprintf('DOF \t err \t\t order \n')
% 
% err0 = 0;
% nU0 = 0;
% for k = 1:length(all_Mesh)
%     
%     Mesh = all_Mesh{k};
%     
%     u_exact = zeros(Mesh.nU, 1);
%     for U = 1:Mesh.nU
%         xc = Mesh.xc(U); yc = Mesh.yc(U);
%         u_exact(U) = PDE.u(xc, yc);
%     end
%     
%     weight = order2_weight(Mesh, PDE);
%     [A, F] = NPS(Mesh, PDE, weight);
%     u = A \ F;
%     
%     nU = Mesh.nU;
%     err = max(abs(u - u_exact));
%     fprintf('%d \t %.2e \t %.2f \n', nU, err, ...
%         -2 * (log(err0)-log(err)) / (log(nU0)-log(nU)));
%     nU0 = nU; err0 = err;
% end

%% FPS
% fprintf('\nFPS:\n')
% fprintf('DOF \t err \t Nit \t order \n')
% 
% err0 = 0;
% nU0 = 0;
% for k = 1:length(all_Mesh)
%     Mesh = all_Mesh{k};
%     
%     u_exact = zeros(Mesh.nU, 1);
%     for U = 1:Mesh.nU
%         xc = Mesh.xc(U); yc = Mesh.yc(U);
%         u_exact(U) = PDE.u(xc, yc);
%     end
%     
%     % picard iteration
%     u0 = zeros(Mesh.nU, 1);
%     for nit = 1:300
%         weight = order2_weight(Mesh, PDE, u0);
%         up = interp_by_weight(Mesh, PDE, u0, weight);
%         up = interp_limit(Mesh, up, u0);
%         [A, F] = FPS(Mesh, PDE, up, u0);
%         u1 = A \ F;
% 
%         if norm(u0-u1) < 1e-10
%             break
%         else
%             u0 = u1;
%         end
%     end
%     
%     nU = Mesh.nU;
%     err = max(abs(u1 - u_exact));
%     fprintf('%d \t %.2e \t %d \t %.2f \n', nU, err, nit, ...
%         -2 * (log(err0)-log(err)) / (log(nU0)-log(nU)));
%     nU0 = nU; err0 = err;
% end

%% EBS1
% fprintf('\nECS-I:\n')
% fprintf('DOF \t err \t order \t ppr err \t order\n')
% 
% derr0 = 0;
% err0 = 0;
% nE0 = 0;
% for k = 1:length(all_Mesh)
%     Mesh = all_Mesh{k};
%     
%     u_exact = zeros(Mesh.nE, 1);
%     du_exact = zeros(Mesh.nE, 2);
%     for E = 1:Mesh.nE
%         xc = Mesh.xe(E); yc = Mesh.ye(E);
%         u_exact(E) = PDE.u(xc, yc);
%         du_exact(E,:) = PDE.du(xc, yc);
%     end
%     
%     [A, F] = ECS1(Mesh, PDE);
%     u = A \ F;
%     du = ppr_continous(Mesh, PDE, u);
%     
%     nE = Mesh.nE;
%     err = max(abs(u - u_exact));
%     derr = max(max(abs(du - du_exact)));
%     fprintf('%d \t %.2e \t %.2f \t %.2e \t %.2f\n', ...
%         nE, err, -2 * (log(err0)-log(err)) / (log(nE0)-log(nE)), ...
%         derr, -2 * (log(derr0)-log(derr)) / (log(nE0)-log(nE)) ...
%         );
%     nE0 = nE; err0 = err; derr0 = derr;
% end

%% EBS2
% gamma = 3;
% 
% fprintf('\nECS-II (\\gamma = %g):\n', gamma)
% fprintf('DOF \t err \t order \n')
% 
% derr0 = 0;
% err0 = 0;
% nE0 = 0;
% for k = 1:length(all_Mesh)
%     Mesh = all_Mesh{k};
%     
%     u_exact = zeros(Mesh.nE, 1);
%     du_exact = zeros(Mesh.nE, 2);
%     for E = 1:Mesh.nE
%         xc = Mesh.xe(E); yc = Mesh.ye(E);
%         u_exact(E) = PDE.u(xc, yc);
%         du_exact(E,:) = PDE.du(xc, yc);
%     end
%     
%     [A, F] = ECS2(Mesh, PDE, gamma);
%     u = A \ F;
%     du = ppr_continous(Mesh, PDE, u);
%     
%     nE = Mesh.nE;
%     err = max(abs(u - u_exact));
%     derr = max(max(abs(du - du_exact)));
%     fprintf('%d \t %.2e \t %.2f \t %.2e \t %.2f\n', ...
%         nE, err, -2 * (log(err0)-log(err)) / (log(nE0)-log(nE)), ...
%         derr, -2 * (log(derr0)-log(derr)) / (log(nE0)-log(nE)) ...
%         );
%     nE0 = nE; err0 = err; derr0 = derr;
% end

%% PPR
fprintf('\nPPR:\n')
fprintf('DOF \t err \t order \n')

derr0 = 0;
nE0 = 0;
for k = 1:length(all_Mesh)
    Mesh = all_Mesh{k};
    
    u_exact = zeros(Mesh.nE, 1);
    du_exact = zeros(Mesh.nE, 2);
    for E = 1:Mesh.nE
        xc = Mesh.xe(E); yc = Mesh.ye(E);
        u_exact(E) = PDE.u(xc, yc);
        du_exact(E,:) = PDE.du(xc, yc);
    end
    
    du = ppr_continous(Mesh, PDE, u_exact);
    
    nE = Mesh.nE;
    derr = max(max(abs(du - du_exact)));
    fprintf('%d \t %.2e \t %.2f\n', ...
        nE, derr, -2 * (log(derr0)-log(derr)) / (log(nE0)-log(nE)) ...
        );
    nE0 = nE; derr0 = derr;
end