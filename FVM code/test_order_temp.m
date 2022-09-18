clear;

old_path = path;
old_path = path(old_path, './mesh');
old_path = path(old_path, './PDE');
old_path = path(old_path, './interp');
old_path = path(old_path, './scheme');
old_path = path(old_path, './utils');

PDE = linear_stab1();

all_Mesh = cell(6,1);
all_Mesh{1} = load_mesh_file('mesh/mesh1_1.txt');
all_Mesh{2} = load_mesh_file('mesh/mesh1_2.txt');
all_Mesh{3} = load_mesh_file('mesh/mesh1_3.txt');
all_Mesh{4} = load_mesh_file('mesh/mesh1_4.txt');
all_Mesh{5} = load_mesh_file('mesh/mesh1_5.txt');
all_Mesh{6} = load_mesh_file('mesh/mesh1_6.txt');

% p = 0.4;
% all_Mesh = cell(6,1);
% all_Mesh{1} = get_Kershaw_mesh(4, 4, p);
% all_Mesh{2} = get_Kershaw_mesh(8, 8, p);
% all_Mesh{3} = get_Kershaw_mesh(16, 16, p);
% all_Mesh{4} = get_Kershaw_mesh(32, 32, p);
% all_Mesh{5} = get_Kershaw_mesh(64, 64, p);
% all_Mesh{6} = get_Kershaw_mesh(128, 128, p);

% p = 0.4;
% all_Mesh = cell(6,1);
% all_Mesh{1} = get_sin_mesh(4, 4, p);
% all_Mesh{2} = get_sin_mesh(8, 8, p);
% all_Mesh{3} = get_sin_mesh(16, 16, p);
% all_Mesh{4} = get_sin_mesh(32, 32, p);
% all_Mesh{5} = get_sin_mesh(64, 64, p);
% all_Mesh{6} = get_sin_mesh(128, 128, p);

gamma = 1;

nE0 = 0;
de00 = 0;
e10 = 0;
de10 = 0;
e20 = 0;
de20 = 0;
for k = 1:length(all_Mesh)
    Mesh = all_Mesh{k};
    
    u_exact = zeros(Mesh.nE, 1);
    for E = 1:Mesh.nE
        xc = Mesh.xe(E); yc = Mesh.ye(E);
        u_exact(E) = PDE.u(xc, yc);
    end
    
    du_exact = zeros(Mesh.nP, 2);
    for P = 1:Mesh.nP
        xc = Mesh.xp(P); yc = Mesh.yp(P);
        du_exact(P,:) = PDE.du(xc, yc);
    end
    
    du0 = ppr_point_temp(Mesh, PDE, u_exact);
    
    [A, F] = ECS1(Mesh, PDE);
    u1 = A \ F;
    du1 = ppr_point_temp(Mesh, PDE, u1);
    
    [A, F] = ECS2(Mesh, PDE, gamma);
    u2 = A \ F;
    du2 = ppr_point_temp(Mesh, PDE, u2);
    
    nE = Mesh.nE;
    de0 = max(max(abs(du0 - du_exact)));
    e1 = max(abs(u1 - u_exact));
    de1 = max(max(abs(du1 - du_exact)));
    e2 = max(abs(u2 - u_exact));
    de2 = max(max(abs(du2 - du_exact)));
    
    fprintf('%d & %.2e & %.2f & %.2e & %.2f & %.2e & %.2f \\\\\n', nE, ...
        de0, -2 * (log(de00)-log(de0)) / (log(nE0)-log(nE)), ...
        e1, -2 * (log(e10)-log(e1)) / (log(nE0)-log(nE)), ...
        de1, -2 * (log(de10)-log(de1)) / (log(nE0)-log(nE)), ...
        e2, -2 * (log(e20)-log(e2)) / (log(nE0)-log(nE)), ...
        de2, -2 * (log(de20)-log(de2)) / (log(nE0)-log(nE)) ...
        );
    
    nE0 = nE;
    de00 = de0;
    e10 = e1;
    de10 = de1;
    e20 = e2;
    de20 = de2;
end
