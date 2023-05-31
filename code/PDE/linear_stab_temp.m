function PDE = linear_stab_temp(delta)

if nargin < 1
    delta = 1e3;
end

    function u = u(x, y)
        u = sin(pi*x) * sin(pi*y);
    end

    function du = du(x, y)
        du = pi * [cos(pi*x)*sin(pi*y); cos(pi*y)*sin(pi*x)];
    end

    function a = a(~, ~, ~)
        a = [delta, 0; 0, 1];
    end

    function f = f(x, y, ~)
        f = (delta+1)*pi*pi * sin(pi*x) * sin(pi*y);
    end

PDE = struct('a', @a, 'u', @u, 'f', @f, 'du', @du, ...
    'bdtype', 'D');
end