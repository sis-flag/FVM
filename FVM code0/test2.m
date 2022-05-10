err0 = 0;
nU0 = 0;
fprintf('N \t err \t order \n')
for N = [8, 16, 32, 64, 128]
    
    Mesh = get_sin_mesh(N, N, 0.5);
    
    PDE = problem5();
    
    weight = interp_order2_weight(Mesh, PDE);
    
    uc = solve_NPS(Mesh, PDE, weight);
    
    err = get_L2_err(Mesh, PDE, uc);
    
    nU = Mesh.nU;
    fprintf('%d \t %.2e \t %g \n', N, err, ...
        -2 * (log(err0)-log(err)) / (log(nU0)-log(nU)));
    nU0 = nU; err0 = err;
end