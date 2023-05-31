clear;
addpath("mesh", "PDE", "scheme", "interp", "utils")


PDE = linear_stab1();

all_Mesh = cell(6,1);
all_Mesh{1} = load_mesh_file('mesh/mesh4_1_1.txt');
all_Mesh{2} = load_mesh_file('mesh/mesh4_1_2.txt');
all_Mesh{3} = load_mesh_file('mesh/mesh4_1_3.txt');
all_Mesh{4} = load_mesh_file('mesh/mesh4_1_4.txt');
all_Mesh{5} = load_mesh_file('mesh/mesh4_1_5.txt');
all_Mesh{6} = load_mesh_file('mesh/mesh4_1_6.txt');

nU = zeros(length(all_Mesh), 1, 'double');
nE = zeros(length(all_Mesh), 1, 'double');
for k = 1:length(all_Mesh)
    nU(k) = all_Mesh{k}.nU;
    nE(k) = all_Mesh{k}.nE;
end

%% NPS
NPS_Linf = zeros(length(all_Mesh), 1, 'double');
NPS_L2 = zeros(length(all_Mesh), 1, 'double');
NPS_time = zeros(length(all_Mesh), 1, 'double');
NPS_nnz = zeros(length(all_Mesh), 1, 'double');
for k = 1:length(all_Mesh)
    Mesh = all_Mesh{k};
    
    u_exact = zeros(Mesh.nU, 1);
    for U = 1:Mesh.nU
        xc = Mesh.xc(U); yc = Mesh.yc(U);
        u_exact(U) = PDE.u(xc, yc);
    end
    
    weight = geo_order2_weight(Mesh, PDE);
    [A, F] = NPS(Mesh, PDE, weight);

    tic;
    u = A \ F;
    NPS_time(k) = toc;

    NPS_nnz(k) = nnz(A);
    NPS_Linf(k) = norm_unit(Mesh, u - u_exact, inf) / ...
        norm_unit(Mesh, u_exact, inf);
    NPS_L2(k) = norm_unit(Mesh, u - u_exact, 2) / ...
        norm_unit(Mesh, u_exact, 2);
end

%% ECS1
ECS1_Linf = zeros(length(all_Mesh), 1, 'double');
ECS1_L2 = zeros(length(all_Mesh), 1, 'double');
ECS1_time = zeros(length(all_Mesh), 1, 'double');
ECS1_nnz = zeros(length(all_Mesh), 1, 'double');
for k = 1:length(all_Mesh)
    Mesh = all_Mesh{k};
    
    u_exact = zeros(Mesh.nE, 1);
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
end

%% ECS2
ECS2_Linf = zeros(length(all_Mesh), 1, 'double');
ECS2_L2 = zeros(length(all_Mesh), 1, 'double');
ECS2_time = zeros(length(all_Mesh), 1, 'double');
ECS2_nnz = zeros(length(all_Mesh), 1, 'double');
for k = 1:length(all_Mesh)
    Mesh = all_Mesh{k};
    
    u_exact = zeros(Mesh.nE, 1);
    for E = 1:Mesh.nE
        xc = Mesh.xe(E); yc = Mesh.ye(E);
        u_exact(E) = PDE.u(xc, yc);
    end
    
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

%% plot error
figure
hold on
plot(nU, NPS_L2, '-^', ...
    'Color', [0 0.4470 0.7410], 'MarkerFaceColor', [0 0.4470 0.7410])
plot(nE, ECS1_L2, '-o', ...
    'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor', [0.8500 0.3250 0.0980])
plot(nE, ECS2_L2, '-s', ...
    'Color', [0.4940 0.1840 0.5560], 'MarkerFaceColor', [0.4940 0.1840 0.5560])
plot(nU, NPS_Linf, '-^', ...
    'Color', [0 0.4470 0.7410], 'MarkerFaceColor', 'none')
plot(nE, ECS1_Linf, '-o', ...
    'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor', 'none')
plot(nE, ECS2_Linf, '-s', ...
    'Color', [0.4940 0.1840 0.5560], 'MarkerFaceColor', 'none')
plot([5e2, 2e4], 2e1*[5e2, 2e4].^(-1), 'k--')
xlabel("nkw")
ylabel("error")
xlim([2e2, 3e4])
ylim([3e-4, 1e-1])
set(gca, 'yscale', 'log')
set(gca, 'xscale', 'log')
grid on
legend(["L^2 error of NPS", "L^2 error of ECS-FV", "L^2 error of ECS-MFD", ...
        "L^\infty error of NPS", "L^\infty error of ECS-FV", "L^\infty error of ECS-MFD", ...
        "order = 2"], ...
        'Location', 'southwest')
set(gcf, 'Position', [200, 200, 400, 400])
print('pics/e3.eps', '-depsc', '-loose')

%% rate
% NPS_L2_Rate = zeros(length(all_Mesh), 1, 'double');
% NPS_Linf_Rate = zeros(length(all_Mesh), 1, 'double');
% ECS1_L2_Rate = zeros(length(all_Mesh), 1, 'double');
% ECS1_Linf_Rate = zeros(length(all_Mesh), 1, 'double');
% ECS2_L2_Rate = zeros(length(all_Mesh), 1, 'double');
% ECS2_Linf_Rate = zeros(length(all_Mesh), 1, 'double');
% for k = 1:length(all_Mesh)
%     if k == 1
%         NPS_L2_Rate(k) = NaN;
%         NPS_Linf_Rate(k) = NaN;
%         ECS1_L2_Rate(k) = NaN;
%         ECS1_Linf_Rate(k) = NaN;
%         ECS2_L2_Rate(k) = NaN;
%         ECS2_Linf_Rate(k) = NaN;
%     else
%         NPS_L2_Rate(k) = -2 * (log(NPS_L2(k-1))-log(NPS_L2(k))) / ...
%             (log(all_Mesh{k-1}.nU)-log(all_Mesh{k}.nU));
%         NPS_Linf_Rate(k) = -2 * (log(NPS_Linf(k-1))-log(NPS_Linf(k))) /...
%             (log(all_Mesh{k-1}.nU)-log(all_Mesh{k}.nU));
%         ECS1_L2_Rate(k) = -2 * (log(ECS1_L2(k-1))-log(ECS1_L2(k))) / ...
%             (log(all_Mesh{k-1}.nE)-log(all_Mesh{k}.nE));
%         ECS1_Linf_Rate(k) = -2 * (log(ECS1_Linf(k-1))-log(ECS1_Linf(k))) /...
%             (log(all_Mesh{k-1}.nE)-log(all_Mesh{k}.nE));
%         ECS2_L2_Rate(k) = -2 * (log(ECS2_L2(k-1))-log(ECS2_L2(k))) / ...
%             (log(all_Mesh{k-1}.nE)-log(all_Mesh{k}.nE));
%         ECS2_Linf_Rate(k) = -2 * (log(ECS2_Linf(k-1))-log(ECS2_Linf(k))) /...
%             (log(all_Mesh{k-1}.nE)-log(all_Mesh{k}.nE));
%     end
% end

%% table
% disp('L2 error')
% for k = 1:length(all_Mesh)
%     if k == 1
%         fprintf('$%d \\times %d$ & %d & %.1e & %s & %d & %.1e & %s & %d & %.1e & %s \\\\\n', ...
%             sqrt(all_Mesh{k}.nU), sqrt(all_Mesh{k}.nU),...
%             all_Mesh{k}.nU, NPS_L2(k), '*', ...
%             all_Mesh{k}.nE, ECS1_L2(k), '*', ...
%             all_Mesh{k}.nE, ECS2_L2(k), '*')
%     else
%         fprintf('$%d \\times %d$ & %d & %.1e & %.2f & %d & %.1e & %.2f & %d & %.1e & %.2f \\\\\n', ...
%             sqrt(all_Mesh{k}.nU), sqrt(all_Mesh{k}.nU),...
%             all_Mesh{k}.nU, NPS_L2(k), NPS_L2_Rate(k), ...
%             all_Mesh{k}.nE, ECS1_L2(k), ECS1_L2_Rate(k), ...
%             all_Mesh{k}.nE, ECS2_L2(k), ECS2_L2_Rate(k))
%     end
% end
% 
% disp('Linf error')
% for k = 1:length(all_Mesh)
%     if k == 1
%         fprintf('$%d \\times %d$ & %d & %.1e & %s & %d & %.1e & %s & %d & %.1e & %s \\\\\n', ...
%             sqrt(all_Mesh{k}.nU), sqrt(all_Mesh{k}.nU),...
%             all_Mesh{k}.nU, NPS_Linf(k), '*', ...
%             all_Mesh{k}.nE, ECS1_Linf(k), '*', ...
%             all_Mesh{k}.nE, ECS2_Linf(k), '*')
%     else
%         fprintf('$%d \\times %d$ & %d & %.1e & %.2f & %d & %.1e & %.2f & %d & %.1e & %.2f \\\\\n', ...
%             sqrt(all_Mesh{k}.nU), sqrt(all_Mesh{k}.nU),...
%             all_Mesh{k}.nU, NPS_Linf(k), NPS_Linf_Rate(k), ...
%             all_Mesh{k}.nE, ECS1_Linf(k), ECS1_Linf_Rate(k), ...
%             all_Mesh{k}.nE, ECS2_Linf(k), ECS2_Linf_Rate(k))
%     end
% end

for k = 1:length(all_Mesh)
fprintf("$%d \\times %d$ & ", sqrt(nU(k)), sqrt(nU(k)))
fprintf("%d & %d & %.1f ms & ", nU(k), NPS_nnz(k), NPS_time(k)*1e3)
fprintf("%d & %d & %.1f ms & ", nE(k), ECS1_nnz(k), ECS1_time(k)*1e3)
fprintf("%d & %d & %.1f ms \\\\\n", nE(k), ECS2_nnz(k), ECS2_time(k)*1e3)
end
