function Mesh = arrange_quad_tri(x, y)
% generate a triangle mesh struct from quadrilateral mesh
% input:
%     X, Y (Nx x Ny array): coordinates of each points.
% output: a mesh struct

[Nx, Ny] = size(x);
Nx = Nx-1; Ny = Ny-1;
coord = [x(:)'; y(:)'];

U2P = cell(2*Nx*Ny, 1);
for ny = 1:Ny
    for nx = 1:Nx
        k = (ny-1)*Nx + nx;
        U2P{2*k-1} = [k+ny+Nx, k+ny-1, k+ny];
        U2P{2*k} = [k+ny+Nx k+ny, k+ny+Nx+1];
    end
end

Mesh = arrange_polygonal(coord, U2P);

end