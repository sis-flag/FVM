clear;

old_path = path;
old_path = path(old_path, './mesh');
old_path = path(old_path, './PDE');
old_path = path(old_path, './interp');

PDE = linear_stab1();

% all_Mesh = cell(5,1);
% all_Mesh{1} = load_mesh_file('mesh/mesh1_1.txt');
% all_Mesh{2} = load_mesh_file('mesh/mesh1_2.txt');
% all_Mesh{3} = load_mesh_file('mesh/mesh1_3.txt');
% all_Mesh{4} = load_mesh_file('mesh/mesh1_4.txt');
% all_Mesh{5} = load_mesh_file('mesh/mesh1_5.txt');

p = 1;
all_Mesh = cell(5,1);
all_Mesh{1} = get_rand_mesh(5, 5, p);
all_Mesh{2} = get_rand_mesh(10, 10, p);
all_Mesh{3} = get_rand_mesh(20, 20, p);
all_Mesh{4} = get_rand_mesh(40, 40, p);
all_Mesh{5} = get_rand_mesh(80, 80, p);

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
% all_Mesh{1} = get_Kershaw_mesh(4, 4, 0.25);
% all_Mesh{2} = get_Kershaw_mesh(8, 8, 0.25);
% all_Mesh{3} = get_Kershaw_mesh(16, 16, 0.25);
% all_Mesh{4} = get_Kershaw_mesh(32, 32, 0.25);
% all_Mesh{5} = get_Kershaw_mesh(64, 64, 0.25);
% all_Mesh{6} = get_Kershaw_mesh(128, 128, 0.25);

%% NPS
fprintf('\nNPS:\n')
fprintf('DOF & err & order \\\\\n')

err0 = 0;
nU0 = 0;
for k = 1:length(all_Mesh)
    
    Mesh = all_Mesh{k};
    
    u_exact = zeros(Mesh.nU, 1);
    for U = 1:Mesh.nU
        xc = Mesh.xc(U); yc = Mesh.yc(U);
        u_exact(U) = PDE.u(xc, yc);
    end
    
    weight = interp_order2_weight(Mesh, PDE);
    [A, F] = mat_NPS(Mesh, PDE, weight);
    u = A \ F;
    
    err = max(abs(u - u_exact));
    
    nU = Mesh.nU;
    fprintf('%d & %.2e & %g \n', nU, err, ...
        -2 * (log(err0)-log(err)) / (log(nU0)-log(nU)));
    nU0 = nU; err0 = err;
end

%% FPS
% fprintf('\nFPS:\n')
% fprintf('DOF & err && Nit & order \n')
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
%     u0 = zeros(Mesh.nU, 1);
%     for nit = 1:300
%         up = interp_order2_limit(Mesh, PDE, u0);
%         [A, F] = mat_FPS(Mesh, PDE, up, u0);
%         u1 = A \ F;
%         
%         if norm(u0 - u1) < 1e-10
%             break
%         else
%             u0 = u1;
%         end
%     end
%     
%     err = max(abs(u1 - u_exact));
%     
%     nU = Mesh.nU;
%     fprintf('%d & %.2e & %d & %g \n', nU, err, nit, ...
%         -2 * (log(err0)-log(err)) / (log(nU0)-log(nU)));
%     nU0 = nU; err0 = err;
% end

%% EBS
fprintf('\nEBS:\n')
fprintf('DOF & err & order \\\\\n')

err0 = 0;
nE0 = 0;
for k = 1:length(all_Mesh)
    
    Mesh = all_Mesh{k};
    
    u_exact = zeros(Mesh.nE, 1);
    for E = 1:Mesh.nE
        xc = Mesh.xe(E); yc = Mesh.ye(E);
        u_exact(E) = PDE.u(xc, yc);
    end
    
    [A, F] = mat_EBS(Mesh, PDE);
    u = A \ F;
    
    err = max(abs(u - u_exact));
    
    nE = Mesh.nE;
    fprintf('%d & %.2e & %g \n', nE, err, ...
        -2 * (log(err0)-log(err)) / (log(nE0)-log(nE)));
    nE0 = nE; err0 = err;
end
