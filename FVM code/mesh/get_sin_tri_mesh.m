function Mesh = get_sin_tri_mesh(Nx, Ny, p)

% generate uniform mesh
xx = linspace(0, 1, Nx+1);
yy = linspace(0, 1, Ny+1);
[y, x] = meshgrid(yy, xx);

% perturbation
px = p * sin(2*pi*x).*sin(2*pi*y) / (2*pi);
x = x + px;
y = y + px;

Mesh = arrange_quad_tri(x, y);

end