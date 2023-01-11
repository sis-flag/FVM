function Mesh = get_tixing_mesh(Nx, Ny, p)
% generate a tixing mesh with size N on domain [0,1] x [0,1]
% author: sis-flag
% input:
%     Nx, Ny (integer): size of the mesh
%     p (nummber): perturbation of the mesh (0 < p < 1)
% output: a mesh struct

% generate uniform mesh
xx = linspace(0, 1, Nx+1);
yy = linspace(0, 1, Ny+1);
[y, x] = meshgrid(yy, xx);

hx = 1/Nx;
hy = 1/Ny;

% perturbation
y(1:2:end, 2:2:end) = y(1:2:end, 2:2:end) - p * hy;
y(2:2:end, 2:2:end) = y(2:2:end, 2:2:end) + p * hy;

% x(2:2:end, 1:2:end) = x(2:2:end, 1:2:end) - p * hx;
% x(2:2:end, 2:2:end) = x(2:2:end, 2:2:end) + p * hx;

Mesh = arrange_quadrilateral(x, y);

end