function plot_func_unit(Mesh, uc)

hold on
cmap = colormap('jet');

for U = 1:Mesh.nU
    nPs = Mesh.U2P{U};
    xps = zeros(2, length(nPs));
    
    for k = 1:length(nPs)
        xps(:,k) = Mesh.xp{nPs(k)};
    end
    fill3(xps(1,:), xps(2,:), uc(U) * ones(1, length(nPs)), cmap, ...
        'LineStyle', 'None')
end

end