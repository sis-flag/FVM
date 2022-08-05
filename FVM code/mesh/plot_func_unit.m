function plot_func_unit(Mesh, uc)

hold on

cmap = colormap('jet');
temp = (uc - min(uc)) / range(uc) * (size(cmap,1)-1);
color = cmap(ceil(temp + 0.5), :);

for U = 1:Mesh.nU
    nPs = Mesh.U2P{U};
    xps = Mesh.xp(nPs);
    yps = Mesh.yp(nPs);
    
    fill3(xps, yps, uc(U) * ones(1, length(nPs)), ...
        color(U,:), 'LineStyle', 'None')
end

colorbar

end