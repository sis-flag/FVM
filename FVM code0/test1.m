Mesh = load_mesh_file('mesh5.txt');
% plot_mesh(Mesh);

PDE = problem5();

weight = interp_order2_weight(Mesh, PDE);

uc = solve_NPS(Mesh, PDE, weight);

err = get_L2_err(Mesh, PDE, uc);
disp('L2 error')
disp(err)

% uce = zeros(Mesh.nU, 1);
% for U = 1:Mesh.nU
%     xc = Mesh.xc{U};
%     uce(U) = PDE.u(xc(1), xc(2));
% end

figure
plot_func_unit(Mesh, uc)
