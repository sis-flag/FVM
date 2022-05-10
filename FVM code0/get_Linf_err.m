function err = get_Linf_err(Mesh, PDE, uc)

err = 0;
for U = 1:Mesh.nU
    xcu = Mesh.xc{U};
    uce = PDE.u(xcu(1), xcu(2));
    err = max(abs(uce - uc(U)), err);
end

end