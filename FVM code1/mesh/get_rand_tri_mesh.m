function Mesh = get_rand_tri_mesh(Nx, Ny, p)
% generate a random triangle mesh with size N on domain [0,1] x [0,1]
% input:
%     Nx, Ny (integer): size of the mesh
%     p (nummber): perturbation of the mesh (0 < p < 1)
%     when p = 0.0, mesh will reduce to the uniform mesh
% output: a mesh struct

% generate uniform mesh
xx = linspace(0, 1, Nx+1);
yy = linspace(0, 1, Ny+1);
[y, x] = meshgrid(yy, xx);

% perturbation
px = p * (rand(Nx-1, Ny-1)*2-1) / (Nx*2);
py = p * (rand(Nx-1, Ny-1)*2-1) / (Ny*2);
x(2:Nx,2:Ny) = x(2:Nx,2:Ny) + px;
y(2:Nx,2:Ny) = y(2:Nx,2:Ny) + py;

Mesh = arrange_quad_tri(x, y);

end