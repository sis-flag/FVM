function F = getF(f, a, x, u)

N = length(x)-1;
h = 1 / N;
F = h*h*arrayfun(f, x(2:N), u(2:N));

F(1) = F(1) + a((x(1)+x(2))/2, (u(1)+u(2))/2) * u(1);
F(N-1) = F(N-1) + a((x(N)+x(N+1))/2, (u(N)+u(N+1))/2) * u(N+1);

end