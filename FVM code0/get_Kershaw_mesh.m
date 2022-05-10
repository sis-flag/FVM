function Mesh = get_Kershaw_mesh(Nx, Ny, p)
% generate a Kershaw z-mesh with size N on domain [0,1] x [0,1]
% input:
%     Nx, Ny (integer): size of the mesh
%     p (nummber): perturbation of the mesh (0 < p < 1)
%     when p = 0.5, mesh will reduce to the uniform mesh
% output: a mesh struct

x = zeros(Nx+1, Ny+1);
y = zeros(Nx+1, Ny+1);

% divide interval into 5 parts
mm1 = floor(Nx/8) + 1;
mm2 = floor(Nx/4) + 1;
mm4 = Nx + 2 - mm1;
mm3 = mm4 - mm2 + mm1;

% set different mesh height on each part
hy1 = zeros(Nx, 1);
for k = 1: Nx+1
    if k <= mm1
        hy1(k) = p;
    elseif k > mm1 && k <= mm2
        hy1(k) = p+(1-2*p)*(k-mm1)/(mm2-mm1);
    elseif k > mm2 && k <= mm3
        hy1(k) = 1-p+(2*p-1)*(k-mm2)/(mm3-mm2);
    elseif k > mm3 && k <= mm4
        hy1(k) = p+(1-2*p)*(k-mm3)/(mm4-mm3);
    elseif k > mm4
        hy1(k) = 1-p;
    end
end

hy2 = 1 - hy1;
hy1 = hy1 / floor(Ny/2);
if mod(Ny, 2) == 0
    hy2 = hy2 / floor(Ny/2);
else
    hy2 = hy2 / floor((Ny+1)/2);
end

% generate mesh
for k = 1: Ny+1
    x(:,k) = linspace(0,1,Nx+1);
    if k <= Ny/2
        y(:,k) = (k-1)*hy1;
    else
        y(:,k) = 1-(Ny+1-k)*hy2;
    end
end

Mesh = arrange_quadrilateral(x, y);

end