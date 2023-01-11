function Mesh = refine_quad_mesh(Mesh0)

N0 = round(sqrt(Mesh0.nP));

if N0*N0 ~= Mesh0.nP
    error('must  be square!')
end

x0 = reshape(Mesh0.xp, N0, N0);
y0 = reshape(Mesh0.yp, N0, N0);

N = 2 * N0 - 1;
x = zeros(N, N);
y = zeros(N, N);

x(1:2:N, 1:2:N) = x0;
y(1:2:N, 1:2:N) = y0;

x(2:2:N, 1:2:N) = (x0(1:N0-1,:) + x0(2:N0,:))/2;
y(2:2:N, 1:2:N) = (y0(1:N0-1,:) + y0(2:N0,:))/2;

x(1:2:N, 2:2:N) = (x0(:,1:N0-1) + x0(:,2:N0))/2;
y(1:2:N, 2:2:N) = (y0(:,1:N0-1) + y0(:,2:N0))/2;

x(2:2:N, 2:2:N) = (x0(1:N0-1,1:N0-1) + x0(2:N0,1:N0-1) + x0(1:N0-1,2:N0) + x0(2:N0,2:N0))/4;
y(2:2:N, 2:2:N) = (y0(1:N0-1,1:N0-1) + y0(2:N0,1:N0-1) + y0(1:N0-1,2:N0) + y0(2:N0,2:N0))/4;

Mesh = arrange_quadrilateral(x, y);

end