syms x e k b u0
pi = sym(pi);

u = sin(pi*x) + b;
a0 = k * u0 + e;
a = subs(a0, u0, u);

disp(simplify(expand(-diff(a*diff(u,x),x))))

f0 = pi*pi*((2*k*u0-k*b+e)*sin(pi*x) - k);
f = subs(f0, u0, u);

disp(simplify(f == -diff(a*diff(u,x),x)))

disp(u)
disp(a0)
disp(diff(a0, u0))
disp(f0)
disp(diff(f0, u0))