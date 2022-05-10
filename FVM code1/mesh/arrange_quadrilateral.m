function Mesh = arrange_quadrilateral(x, y)
% generate a mesh struct for quadrilateral mesh
% input:
%     X, Y (Nx x Ny array): coordinates of each points.
% output: a mesh struct

[Nx, Ny] = size(x);
Nx = Nx-1; Ny = Ny-1;
coord = [x(:)'; y(:)'];

U2P = cell(Nx*Ny, 1);
for ny = 1:Ny
    for nx = 1:Nx
        k = (ny-1)*Nx + nx;
        U2P{k} = [k+ny+Nx, k+ny-1, k+ny, k+ny+Nx+1];
    end
end

Mesh = arrange_polygonal(coord, U2P);

end
