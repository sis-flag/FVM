clear;

old_path = path;
old_path = path(old_path, './mesh');
old_path = path(old_path, './PDE');
old_path = path(old_path, './interp');
old_path = path(old_path, './scheme');
old_path = path(old_path, './utils');

%% mesh
all_Mesh = cell(6,1);
all_Mesh{1} = load_mesh_file('mesh/mesh3_2.txt');
load mesh/mesh11_3
all_Mesh{2} = arrange_polygonal(coord, U2P);
all_Mesh{3} = load_mesh_file('mesh/mesh4_1_1.txt');
all_Mesh{4} = get_Kershaw_mesh(16, 16, 0.2);
all_Mesh{5} = load_mesh_file('mesh/mesh1_2.txt');
all_Mesh{6} = get_sin_mesh(16, 16, 0.7);

for k = 1:length(all_Mesh)
    figure
    plot_mesh(all_Mesh{k});
    set(gcf, 'Position', [200, 200, 300, 300])
    print(['pics/m', int2str(k),'.eps'],'-depsc','-loose')
end

%% test
PDE = linear_stab0();

ECS1_Linf = zeros(length(all_Mesh), 1, 'double');
ECS1_L2 = zeros(length(all_Mesh), 1, 'double');
ECS2_Linf = zeros(length(all_Mesh), 1, 'double');
ECS2_L2 = zeros(length(all_Mesh), 1, 'double');
for k = 1:length(all_Mesh)
    Mesh = all_Mesh{k};
    
    u_exact = zeros(Mesh.nU, 1);
    for E = 1:Mesh.nE
        xc = Mesh.xe(E); yc = Mesh.ye(E);
        u_exact(E) = PDE.u(xc, yc);
    end

    [A, F] = ECS1(Mesh, PDE);
    u = A \ F;

    ECS1_Linf(k) = norm_edge(Mesh, u - u_exact, inf) / ...
        norm_edge(Mesh, u_exact, inf);
    ECS1_L2(k) = norm_edge(Mesh, u - u_exact, 2) / ...
        norm_edge(Mesh, u_exact, 2);
    
    [A, F] = ECS2(Mesh, PDE);
    u = A \ F;

    ECS2_Linf(k) = norm_edge(Mesh, u - u_exact, inf) / ...
        norm_edge(Mesh, u_exact, inf);
    ECS2_L2(k) = norm_edge(Mesh, u - u_exact, 2) / ...
        norm_edge(Mesh, u_exact, 2);
end

%% table
Tab{1} = '';
Tab{2} = '$L^2$ (ECS-I)';
Tab{3} = '$L^\\infty$ (ECS-I)';
Tab{4} = '$L^2$ (ECS-II)';
Tab{5} = '$L^\\infty$ (ECS-II)';
for k = 1:length(all_Mesh)
    Tab{1} = [Tab{1}, sprintf(' & mesh %d', k)];
    Tab{2} = [Tab{2}, sprintf(' & %.2e', ECS1_L2(k))];
    Tab{3} = [Tab{3}, sprintf(' & %.2e', ECS1_Linf(k))];
    Tab{4} = [Tab{4}, sprintf(' & %.2e', ECS2_L2(k))];
    Tab{5} = [Tab{5}, sprintf(' & %.2e', ECS2_Linf(k))];
end
for k = 1:5
    fprintf([Tab{k}, ' \\\\\n'])
end