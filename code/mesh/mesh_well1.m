Nx = 54; Ny = 54;

% mesh points
xx = linspace(0, 1, Nx+1);
yy = linspace(0, 1, Ny+1);
[y, x] = meshgrid(yy, xx);

% perturbation
px = (rand(Nx-1, Ny-1)*2-1) / 180;
py = (rand(Nx-1, Ny-1)*2-1) / 180;
px([24,30], 24:30) = 0;
py([24,30], 24:30) = 0;
px(24:30, [24,30]) = 0;
py(24:30, [24,30]) = 0;
x(2:Nx,2:Ny) = x(2:Nx,2:Ny) + px;
y(2:Nx,2:Ny) = y(2:Nx,2:Ny) + py;

% delete hole
x(26:30, 26:30) = nan;
y(26:30, 26:30) = nan;

x = x(:);
y = y(:);

% indxes in hole
ind = find(~isnan(x));
fuck = nan * ones(length(x), 1);
for k = ind
    fuck(ind) = k;
end


k = 0;
U2P = cell(Nx*Ny-6*6, 1);
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
save_mesh_file("mesh_well1.txt", Mesh)