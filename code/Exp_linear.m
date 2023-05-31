clear;
addpath("mesh", "PDE", "scheme", "interp", "utils")
rng(0);

%% mesh
all_Mesh = cell(8,1);
all_Mesh{1} = load_mesh_file('mesh/mesh3_2.txt');
load mesh/mesh11_3
all_Mesh{2} = arrange_polygonal(coord, U2P);
all_Mesh{3} = load_mesh_file('mesh/mesh4_1_1.txt');
all_Mesh{4} = get_Kershaw_mesh(16, 16, 0.2);
all_Mesh{5} = load_mesh_file('mesh/mesh1_2.txt');
all_Mesh{6} = get_rand_mesh(16, 16, 0.7);
all_Mesh{7} = load_mesh_file('mesh/mesh7.txt');
all_Mesh{8} = get_sin_mesh(16, 16, 0.7);

% for k = 1:length(all_Mesh)
%     figure
%     plot_mesh(all_Mesh{k});
%     axis("equal")
%     xlim([0, 1])
%     ylim([0, 1])
%     set(gcf, 'Position', [200, 200, 300, 300])
%     print(['pics/m', int2str(k),'.eps'],'-depsc','-loose')
% end

%% problems

all_PDE = cell(8,1);
all_PDE{1} = linear_stab11();
all_PDE{2} = linear_stab0();
all_PDE{3} = linear_stab0();
all_PDE{4} = linear_stab0();
all_PDE{5} = linear_stab11();
all_PDE{6} = linear_stab0();
all_PDE{7} = linear_stab7();
all_PDE{8} = linear_stab11();

%% test
ECS1_Linf = zeros(length(all_Mesh), 1, 'double');
ECS1_L2 = zeros(length(all_Mesh), 1, 'double');
ECS1_time = zeros(length(all_Mesh), 1, 'double');
ECS1_nnz = zeros(length(all_Mesh), 1, 'double');
ECS2_Linf = zeros(length(all_Mesh), 1, 'double');
ECS2_L2 = zeros(length(all_Mesh), 1, 'double');
ECS2_time = zeros(length(all_Mesh), 1, 'double');
ECS2_nnz = zeros(length(all_Mesh), 1, 'double');
for k = 1:length(all_Mesh)
    Mesh = all_Mesh{k};
    PDE = all_PDE{k};
    
    u_exact = zeros(Mesh.nU, 1);
    for E = 1:Mesh.nE
        xc = Mesh.xe(E); yc = Mesh.ye(E);
        u_exact(E) = PDE.u(xc, yc);
    end

    [A, F] = ECS1(Mesh, PDE);

    tic;
    u = A \ F;
    ECS1_time(k) = toc;

    ECS1_nnz(k) = nnz(A);
    ECS1_Linf(k) = norm_edge(Mesh, u - u_exact, inf) / ...
        norm_edge(Mesh, u_exact, inf);
    ECS1_L2(k) = norm_edge(Mesh, u - u_exact, 2) / ...
        norm_edge(Mesh, u_exact, 2);
    
    [A, F] = ECS2(Mesh, PDE);

    tic;
    u = A \ F;
    ECS2_time(k) = toc;

    ECS2_nnz(k) = nnz(A);
    ECS2_Linf(k) = norm_edge(Mesh, u - u_exact, inf) / ...
        norm_edge(Mesh, u_exact, inf);
    ECS2_L2(k) = norm_edge(Mesh, u - u_exact, 2) / ...
        norm_edge(Mesh, u_exact, 2);
end

%% results
for k = 1:length(all_Mesh)
    fprintf("%d & %d & %.1f ms & %.1e & %.1e & %d & %.1f ms & %.1e & %.1e \\\\\n", ...
        all_Mesh{k}.nE, ...
        ECS1_nnz(k), ECS1_time(k)*1e3, ECS1_Linf(k), ECS1_L2(k), ...
        ECS2_nnz(k), ECS2_time(k)*1e3, ECS2_Linf(k), ECS2_L2(k) ...
        )
end