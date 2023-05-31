function pcolor_func_unit(Mesh, uc)

hold on

Mu = max(uc);
mu = min(uc);

cmap = colormap();
ncolor = size(cmap, 1);
temp = (uc - mu) / (Mu - mu) * (ncolor-1) + 1;
color = interp1(1:ncolor, cmap, temp, 'linear');

for U = 1:Mesh.nU
    nPs = Mesh.U2P{U};
    xps = Mesh.xp(nPs);
    yps = Mesh.yp(nPs);
    
    fill(xps, yps, color(U,:), 'LineStyle', '-')
end

colorbar
clim([mu, Mu])

end