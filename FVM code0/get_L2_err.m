function err = get_L2_err(Mesh, PDE, uc)

err = 0;
for U = 1:Mesh.nU
    xcu = Mesh.xc{U};
    uce = PDE.u(xcu(1), xcu(2));
    err = err + Mesh.area{U} * (uce - uc(U))^2;
end
err = sqrt(err);

end