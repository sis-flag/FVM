clear;

old_path = path;
old_path = path(old_path, './mesh');
old_path = path(old_path, './PDE');
old_path = path(old_path, './interp');

PDE = linear_stab6();

p = 0.7;
all_Mesh = cell(9,1);
all_Mesh{1} = get_sin_mesh(8, 8, p);
all_Mesh{2} = get_sin_mesh(11, 11, p);
all_Mesh{3} = get_sin_mesh(16, 16, p);
all_Mesh{4} = get_sin_mesh(23, 23, p);
all_Mesh{5} = get_sin_mesh(32, 32, p);
all_Mesh{6} = get_sin_mesh(45, 45, p);
all_Mesh{7} = get_sin_mesh(64, 64, p);
all_Mesh{8} = get_sin_mesh(90, 90, p);
all_Mesh{9} = get_sin_mesh(128, 128, p);

all_gamma = [0.01, 0.1, 1, 10, 100];

%% test
Linf = zeros(length(all_Mesh), length(all_gamma)+2, 'double');
L2 = zeros(length(all_Mesh), length(all_gamma)+2, 'double');
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
    
    Linf(k, 1) = norm_unit(Mesh, u - u_exact, inf) / ...
        norm_unit(Mesh, u_exact, inf);
    L2(k, 1) = norm_unit(Mesh, u - u_exact, 2) / ...
        norm_unit(Mesh, u_exact, 2);
    
    u_exact = zeros(Mesh.nE, 1);
    for E = 1:Mesh.nE
        xc = Mesh.xe(E); yc = Mesh.ye(E);
        u_exact(E) = PDE.u(xc, yc);
    end
    
    [A, F] = mat_EBS1(Mesh, PDE);
    u = A \ F;
    
    Linf(k, 2) = norm_edge(Mesh, u - u_exact, inf) / ...
        norm_edge(Mesh, u_exact, inf);
    L2(k, 2) = norm_edge(Mesh, u - u_exact, 2) / ...
        norm_edge(Mesh, u_exact, 2);
    
    for kg = 1:length(all_gamma)
        [A, F] = mat_EBS2(Mesh, PDE, all_gamma(kg));
        u = A \ F;
    
        Linf(k,kg+2) = norm_edge(Mesh, u - u_exact, inf) / ...
            norm_edge(Mesh, u_exact, inf);
        L2(k,kg+2) = norm_edge(Mesh, u - u_exact, 2) / ...
            norm_edge(Mesh, u_exact, 2);
    end
end

%% plot
nUs = zeros(length(all_Mesh), 1, 'double');
nEs = zeros(length(all_Mesh), 1, 'double');
for k = 1:length(all_Mesh)
    nUs(k) = all_Mesh{k}.nU;
    nEs(k) = all_Mesh{k}.nE;
end

leg = cell(length(all_gamma)+3, 1);
leg{1} = 'slope=-2';
leg{2} = 'NPS';
leg{3} = 'ECS-I';
for k = 1:length(all_gamma)
    leg{k+3} = sprintf('ECS-II(\\gamma=%g)', all_gamma(k));
end

figure
set(gcf, 'Position', [200, 200, 1000, 400])

subplot(1,2,1)
hold on
plot([8, 180], 7*[8^(-2), 180^(-2)], 'k--')
plot(sqrt(nUs), Linf(:,1), 'm-*');
for k = 2:length(all_gamma)+2
    plot(sqrt(nEs), Linf(:,k), '-*');
end
grid('on')
xlim([5, 3e2])
ylim([8e-5, 8e0])
xlabel('sqrt(nkw)')
ylabel('L^\infty error')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'FontSize', 12)
set(gca, 'Position', [0.1, 0.15, 0.3, 0.8])

subplot(1,2,2)
hold on
plot([8, 180], 3*[8^(-2), 180^(-2)], 'k--')
plot(sqrt(nUs), L2(:,1), 'm-*');
for k = 2:length(all_gamma)+2
    plot(sqrt(nEs), L2(:,k), '-*');
end
grid('on')
xlim([5, 3e2])
ylim([5e-5, 5e0])
xlabel('sqrt(nkw)')
ylabel('L^2 error')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'FontSize', 12)
set(gca, 'Position', [0.5, 0.15, 0.3, 0.8])

o_o = legend(leg);
set(o_o, 'Position', [0.85, 0.5, 0.1, 0.4])
print('pics/m6err.eps','-depsc','-loose')
