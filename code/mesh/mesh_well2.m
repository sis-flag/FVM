Nx = 33; Ny = 33;

% mesh points
xx = linspace(0, 1, Nx+1);
yy = linspace(0, 1, Ny+1);
[y, x] = meshgrid(yy, xx);

% delete hole
x(11:12, 17:18) = nan;
y(11:12, 17:18) = nan;
x(23:24, 17:18) = nan;
y(23:24, 17:18) = nan;

x = x(:);
y = y(:);

% indxes in hole
ind = find(~isnan(x));
fuck = nan * ones(length(x), 1);
for k = ind
    fuck(ind) = k;
end


k = 0;
U2P = cell(Nx*Ny-9*9*2, 1);
for ny = 1:Ny
    for nx = 1:Nx
        temp = [fuck((ny-1)*Nx + nx + ny + Nx), ...
                fuck((ny-1)*Nx + nx + ny - 1), ...
                fuck((ny-1)*Nx + nx + ny), ...
                fuck((ny-1)*Nx + nx + ny  + Nx + 1)
                ];
        if all(temp > 0)
            k = k+1;
            U2P{k} = temp;
        end
    end
end

coord = [x(:)'; y(:)'];
Mesh = arrange_polygonal(coord, U2P);

plot_mesh(Mesh)
save_mesh_file("mesh_well2.txt", Mesh)