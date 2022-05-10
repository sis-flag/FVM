v = [0,8;8,2;5,-8;-5,-8;-8,2];
[x,y] = meshgrid(-15:0.5:15, -15:0.5:15);
phi = zeros([size(v,1), size(x)]);
for i = 1:size(x,1)
    for j = 1:size(x,2)
        phi(:,i,j) = mv2dglob(v, [x(i,j), y(i,j)]);
    end
end

rphi = real(phi);

figure
hold on

o_o = zeros(size(v,1),1);
o_o(1) = 1;
for k = 1:size(v,1)-1
    plot3(v(k:k+1,1), v(k:k+1,2), o_o(k:k+1), 'o-');
end
plot3(v([end,1],1), v([end,1],2), o_o([end,1]), 'o-');

cl = surf(x, y, reshape(rphi(1,:,:), size(x)));
cl.LineStyle = 'None';

