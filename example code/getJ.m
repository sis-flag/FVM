function J = getJ(a, au, fu, x, u)

N = length(x)-1;
h = 1 / N;
xm = conv(x, [0.5, 0.5], 'valid');
um = conv(u, [0.5, 0.5], 'valid');
am = arrayfun(a, xm, um);
aum = arrayfun(au, xm, um);
fuu = h*h*arrayfun(fu, x(2:N), u(2:N));

J = zeros(N-1);
for k = 1:N-1
    J(k,k) = am(k) + am(k+1) - fuu(k) ...
        + 0.5 * aum(k) * (u(k+1) - u(k)) ...
        + 0.5 * aum(k+1) * (u(k+1) - u(k+2));
end
for k = 1:N-2
    J(k,k+1) = -am(k+1) + 0.5 * aum(k+1) * (u(k+1) - u(k+2));
    J(k+1,k) = -am(k+1) + 0.5 * aum(k+1) * (u(k+2) - u(k+1));
end
% 
% 
% for k =1:N-1
%     D(k,k) = am(k)...
%         + a((x(k+1)+x(k+2))/2, (u(k+1)+u(k+2))/2)...
%         + 0.5 * au((x(k+1)+x(k))/2, (u(k+1)+u(k))/2) * (u(k+1) - u(k))...
%         + 0.5 * au((x(k+1)+x(k+2))/2, (u(k+1)+u(k+2))/2) * (u(k+1) - u(k+2));
%     if k < N-1
%         D(k,k+1) = -a((x(k+1)+x(k+2))/2, (u(k+1)+u(k+2))/2)...
%             + 0.5 * au((x(k+1)+x(k+2))/2, (u(k+1)+u(k+2))/2) * (u(k+1) - u(k+2));
%     end
%     if k > 1
%         D(k,k-1) = -a((x(k+1)+x(k))/2, (u(k+1)+u(k))/2)...
%             + 0.5 * au((x(k+1)+x(k))/2, (u(k+1)+u(k))/2) * (u(k+1) - u(k));
%     end
% end

end