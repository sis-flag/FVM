clear;

old_path = path;
old_path = path(old_path, './mesh');
old_path = path(old_path, './PDE');
old_path = path(old_path, './interp');

PDE = linear_stab3();

p = 0.2;
all_Mesh = cell(6,1);
all_Mesh{1} = get_Kershaw_mesh(4, 4, p);
all_Mesh{2} = get_Kershaw_mesh(8, 8, p);
all_Mesh{3} = get_Kershaw_mesh(16, 16, p);
all_Mesh{4} = get_Kershaw_mesh(32, 32, p);
all_Mesh{5} = get_Kershaw_mesh(64, 64, p);
all_Mesh{6} = get_Kershaw_mesh(128, 128, p);

%% NPS
NPS_Linf = zeros(length(all_Mesh), 1, 'double');
NPS_L2 = zeros(length(all_Mesh), 1, 'double');
for k = 1:length(all_Mesh)
    Mesh = all_Mesh{k};
    
    u_exact = zeros(Mesh.nU, 1);
    for U = 1:Mesh.nU
        xc = Mesh.xc(U); yc = Mesh.yc(U);
        u_exact(U) = PDE.u(xc, yc);
    end
    
    weight = order2_weight(Mesh, PDE);
    [A, F] = mat_NPS(Mesh, PDE, weight);
    u = A \ F;
    
    NPS_Linf(k) = norm_unit(Mesh, u - u_exact, inf) / ...
        norm_unit(Mesh, u_exact, inf);
    NPS_L2(k) = norm_unit(Mesh, u - u_exact, 2) / ...
        norm_unit(Mesh, u_exact, 2);
end

%% EBS1
EBS1_Linf = zeros(length(all_Mesh), 1, 'double');
EBS1_L2 = zeros(length(all_Mesh), 1, 'double');
for k = 1:length(all_Mesh)
    Mesh = all_Mesh{k};
    
    u_exact = zeros(Mesh.nE, 1);
    for E = 1:Mesh.nE
        xc = Mesh.xe(E); yc = Mesh.ye(E);
        u_exact(E) = PDE.u(xc, yc);
    end
    
    [A, F] = mat_EBS1(Mesh, PDE);
    u = A \ F;
    
    EBS1_Linf(k) = norm_edge(Mesh, u - u_exact, inf) / ...
        norm_edge(Mesh, u_exact, inf);
    EBS1_L2(k) = norm_edge(Mesh, u - u_exact, 2) / ...
        norm_edge(Mesh, u_exact, 2);
end

%% EBS2
EBS2_Linf = zeros(length(all_Mesh), 1, 'double');
EBS2_L2 = zeros(length(all_Mesh), 1, 'double');
for k = 1:length(all_Mesh)
    Mesh = all_Mesh{k};
    
    u_exact = zeros(Mesh.nE, 1);
    for E = 1:Mesh.nE
        xc = Mesh.xe(E); yc = Mesh.ye(E);
        u_exact(E) = PDE.u(xc, yc);
    end
    
    [A, F] = mat_EBS2(Mesh, PDE);
    u = A \ F;
    
    EBS2_Linf(k) = norm_edge(Mesh, u - u_exact, inf) / ...
        norm_edge(Mesh, u_exact, inf);
    EBS2_L2(k) = norm_edge(Mesh, u - u_exact, 2) / ...
        norm_edge(Mesh, u_exact, 2);
end

%% rate
NPS_L2_Rate = zeros(length(all_Mesh), 1, 'double');
NPS_Linf_Rate = zeros(length(all_Mesh), 1, 'double');
EBS1_L2_Rate = zeros(length(all_Mesh), 1, 'double');
EBS1_Linf_Rate = zeros(length(all_Mesh), 1, 'double');
EBS2_L2_Rate = zeros(length(all_Mesh), 1, 'double');
EBS2_Linf_Rate = zeros(length(all_Mesh), 1, 'double');
for k = 1:length(all_Mesh)
    if k == 1
        NPS_L2_Rate(k) = NaN;
        NPS_Linf_Rate(k) = NaN;
        EBS1_L2_Rate(k) = NaN;
        EBS1_Linf_Rate(k) = NaN;
        EBS2_L2_Rate(k) = NaN;
        EBS2_Linf_Rate(k) = NaN;
    else
        NPS_L2_Rate(k) = -2 * (log(NPS_L2(k-1))-log(NPS_L2(k))) / ...
            (log(all_Mesh{k-1}.nU)-log(all_Mesh{k}.nU));
        NPS_Linf_Rate(k) = -2 * (log(NPS_Linf(k-1))-log(NPS_Linf(k))) /...
            (log(all_Mesh{k-1}.nU)-log(all_Mesh{k}.nU));
        EBS1_L2_Rate(k) = -2 * (log(EBS1_L2(k-1))-log(EBS1_L2(k))) / ...
            (log(all_Mesh{k-1}.nE)-log(all_Mesh{k}.nE));
        EBS1_Linf_Rate(k) = -2 * (log(EBS1_Linf(k-1))-log(EBS1_Linf(k))) /...
            (log(all_Mesh{k-1}.nE)-log(all_Mesh{k}.nE));
        EBS2_L2_Rate(k) = -2 * (log(EBS2_L2(k-1))-log(EBS2_L2(k))) / ...
            (log(all_Mesh{k-1}.nE)-log(all_Mesh{k}.nE));
        EBS2_Linf_Rate(k) = -2 * (log(EBS2_Linf(k-1))-log(EBS2_Linf(k))) /...
            (log(all_Mesh{k-1}.nE)-log(all_Mesh{k}.nE));
    end
end

%% table
disp('L2 error')
for k = 1:length(all_Mesh)
    if k == 1
        fprintf('$%d \\times %d$ & %d & %.2e & %s & %d & %.2e & %s & %d & %.2e & %s \\\\\n', ...
            sqrt(all_Mesh{k}.nU), sqrt(all_Mesh{k}.nU),...
            all_Mesh{k}.nU, NPS_L2(k), '*', ...
            all_Mesh{k}.nE, EBS1_L2(k), '*', ...
            all_Mesh{k}.nE, EBS2_L2(k), '*')
    else
        fprintf('$%d \\times %d$ & %d & %.2e & %.2f & %d & %.2e & %.2f & %d & %.2e & %.2f \\\\\n', ...
            sqrt(all_Mesh{k}.nU), sqrt(all_Mesh{k}.nU),...
            all_Mesh{k}.nU, NPS_L2(k), NPS_L2_Rate(k), ...
            all_Mesh{k}.nE, EBS1_L2(k), EBS1_L2_Rate(k), ...
            all_Mesh{k}.nE, EBS2_L2(k), EBS2_L2_Rate(k))
    end
end

disp('Linf error')
for k = 1:length(all_Mesh)
    if k == 1
        fprintf('$%d \\times %d$ & %d & %.2e & %s & %d & %.2e & %s & %d & %.2e & %s \\\\\n', ...
            sqrt(all_Mesh{k}.nU), sqrt(all_Mesh{k}.nU),...
            all_Mesh{k}.nU, NPS_Linf(k), '*', ...
            all_Mesh{k}.nE, EBS1_Linf(k), '*', ...
            all_Mesh{k}.nE, EBS2_Linf(k), '*')
    else
        fprintf('$%d \\times %d$ & %d & %.2e & %.2f & %d & %.2e & %.2f & %d & %.2e & %.2f \\\\\n', ...
            sqrt(all_Mesh{k}.nU), sqrt(all_Mesh{k}.nU),...
            all_Mesh{k}.nU, NPS_Linf(k), NPS_Linf_Rate(k), ...
            all_Mesh{k}.nE, EBS1_Linf(k), EBS1_Linf_Rate(k), ...
            all_Mesh{k}.nE, EBS2_Linf(k), EBS2_Linf_Rate(k))
    end
end