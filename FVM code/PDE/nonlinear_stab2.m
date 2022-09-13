function PDE = nonlinear_stab2(e, b, k)

if nargin < 3
    k = 1; e = 0.1; b = 0.1;
end

    function u = u(x, y)
        u = sin(pi*x) * sin(pi*y) + b;
    end

    function du = du(x, y)
        du = pi * [cos(pi*x)*sin(pi*y); cos(pi*y)*sin(pi*x)];
    end

    function a = a(~, ~, u)
        a = (k * u + e) * eye(2);
    end

    function f = f(x, y, u)
        f = 4*k*sin(pi*x)^2*sin(pi*y)^2 ...
            - k*(sin(pi*y)^2 + sin(pi*x)^2) ...
            + 2*(k*b+e)* (u-b);
        f = pi*pi*f;
    end

PDE = struct('a', @a, 'u', @u, 'f', @f, 'du', @du);
end