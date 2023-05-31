clear;
addpath("mesh", "PDE", "scheme", "interp", "utils")

%% Mesh
Mesh = load_mesh_file('mesh/mesh3_2.txt');
disp('locally refined mesh')

%% PDE
PDE = linear_stab2();

u_exact = zeros(Mesh.nU, 1);
for E = 1:Mesh.nE
    xc = Mesh.xe(E); yc = Mesh.ye(E);
    u_exact(E) = PDE.u(xc, yc);
end

%% EBS1
[A, F] = ECS1(Mesh, PDE);
u = A \ F;

L2_err = norm_edge(Mesh, u - u_exact, 2) / norm_edge(Mesh, u_exact, 2);
Linf_err = norm_edge(Mesh, u - u_exact, inf) / norm_edge(Mesh, u_exact, inf);
fprintf('ECS-FV: L2 error %.2e, Linf error %.2e\n', L2_err, Linf_err)

%% EBS1 plot
figure
pcolor_func_edge(Mesh, u)
axis equal
xlim([0, 1])
ylim([0, 1])
set(gcf, 'Position', [200, 200, 400, 330])
print('pics/m1s1u.eps','-depsc','-loose')

figure
pcolor_func_edge(Mesh, u - u_exact)
axis equal
xlim([0, 1])
ylim([0, 1])
set(gcf, 'Position', [200, 200, 400, 330])
print('pics/m1s1err.eps','-depsc','-loose')

%% EBS2
[A, F] = ECS2(Mesh, PDE);
u = A \ F;

L2_err = norm_edge(Mesh, u - u_exact, 2) / norm_edge(Mesh, u_exact, 2);
Linf_err = norm_edge(Mesh, u - u_exact, inf) / norm_edge(Mesh, u_exact, inf);
fprintf('ECS-MFD: L2 error %.2e, Linf error %.2e\n', L2_err, Linf_err)

%% EBS2 plot
figure
pcolor_func_edge(Mesh, u)
axis equal
xlim([0, 1])
ylim([0, 1])
set(gcf, 'Position', [200, 200, 400, 330])
print('pics/m1s2u.eps','-depsc','-loose')

figure
pcolor_func_edge(Mesh, u - u_exact)
axis equal
xlim([0, 1])
ylim([0, 1])
set(gcf, 'Position', [200, 200, 400, 330])
print('pics/m1s2err.eps','-depsc','-loose')

%% Mesh
load mesh/mesh11_3
Mesh = arrange_polygonal(coord, U2P);
disp('random polygonal mesh')

%% PDE
u_exact = zeros(Mesh.nU, 1);
for E = 1:Mesh.nE
    xc = Mesh.xe(E); yc = Mesh.ye(E);
    u_exact(E) = PDE.u(xc, yc);
end

%% EBS1
[A, F] = ECS1(Mesh, PDE);
u = A \ F;

L2_err = norm_edge(Mesh, u - u_exact, 2) / norm_edge(Mesh, u_exact, 2);
Linf_err = norm_edge(Mesh, u - u_exact, inf) / norm_edge(Mesh, u_exact, inf);
fprintf('ECS-FV: L2 error %.2e, Linf error %.2e\n', L2_err, Linf_err)

%% EBS1 plot
figure
pcolor_func_edge(Mesh, u)
axis equal
xlim([0, 1])
ylim([0, 1])
set(gcf, 'Position', [200, 200, 400, 330])
print('pics/m2s1u.eps','-depsc','-loose')

figure
pcolor_func_edge(Mesh, u - u_exact)
axis equal
xlim([0, 1])
ylim([0, 1])
set(gcf, 'Position', [200, 200, 400, 330])
print('pics/m2s1err.eps','-depsc','-loose')

%% EBS2
[A, F] = ECS2(Mesh, PDE);
u = A \ F;

L2_err = norm_edge(Mesh, u - u_exact, 2) / norm_edge(Mesh, u_exact, 2);
Linf_err = norm_edge(Mesh, u - u_exact, inf) / norm_edge(Mesh, u_exact, inf);
fprintf('ECS-MFD: L2 error %.2e, Linf error %.2e\n', L2_err, Linf_err)

%% EBS2 plot
figure
pcolor_func_edge(Mesh, u)
axis equal
xlim([0, 1])
ylim([0, 1])
set(gcf, 'Position', [200, 200, 400, 330])
print('pics/m2s2u.eps','-depsc','-loose')

figure
pcolor_func_edge(Mesh, u - u_exact)
axis equal
xlim([0, 1])
ylim([0, 1])
set(gcf, 'Position', [200, 200, 400, 330])
print('pics/m2s2err.eps','-depsc','-loose')
