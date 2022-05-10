function Mesh = get_sin_mesh(Nx, Ny, p)
% generate a sin mesh with size N on domain [0,1] x [0,1]
% author: sis-flag
% input:
%     Nx, Ny (integer): size of the mesh
%     p (nummber): perturbation of the mesh (0 < p < 1)
%     when p = 0, mesh will reduce to the uniform mesh
% output: a mesh struct

% generate uniform mesh
xx = linspace(0, 1, Nx+1);
yy = linspace(0, 1, Ny+1);
[y, x] = meshgrid(yy, xx);

% perturbation
px = p * sin(2*pi*x).*sin(2*pi*y) / (2*pi);
x = x + px;
y = y + px;

Mesh = arrange_quadrilateral(x, y);

end