function [new_x, new_y] = refine_quadrilateral(old_x, old_y)
% refine old mesh and get new mesh
% input:
%     old_x, old_y: array with size (Nx+1, Ny+1)
% output:
%     new_x, new_y: array with size (2*Nx+1, 2*Ny+1)

Nx = size(old_x, 1) - 1;
Ny = size(old_x, 2) - 1;
new_x = zeros(2*Nx+1, 2*Ny+1);
new_y = zeros(2*Nx+1, 2*Ny+1);

% retain old nodes
new_x(1:2:end, 1:2:end) = old_x;
new_y(1:2:end, 1:2:end) = old_y;

% set new nodes on center of edge
new_x(1:2:end, 2:2:end) = conv2(old_x, [1/2, 1/2], 'valid');
new_y(1:2:end, 2:2:end) = conv2(old_y, [1/2, 1/2], 'valid');
new_x(2:2:end, 1:2:end) = conv2(old_x, [1/2; 1/2], 'valid');
new_y(2:2:end, 1:2:end) = conv2(old_y, [1/2; 1/2], 'valid');

% set new nodes on center of volume
new_x(2:2:end, 2:2:end) = conv2(old_x, [1,1;1,1]/4, 'valid');
new_y(2:2:end, 2:2:end) = conv2(old_y, [1,1;1,1]/4, 'valid');

end