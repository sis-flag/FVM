function A = getA(a, x, u)

N = length(x)-1;
xm = conv(x, [0.5, 0.5], 'valid');
um = conv(u, [0.5, 0.5], 'valid');
am = arrayfun(a, xm, um);

A = zeros(N-1);
for k = 1:N-1
    A(k,k) = am(k) + am(k+1);
end
for k = 1:N-2
    A(k,k+1) = -am(k+1);
    A(k+1,k) = -am(k+1);
end

end

