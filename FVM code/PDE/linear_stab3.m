function PDE = linear_stab3(delta)

if nargin < 1
    delta = 1e5;
end

% syms x y delta
% u = sin(2*pi*x) * exp(-2*pi/sqrt(delta)*y)
% a = [1, 0; 0, delta]
% F = a * [diff(u, x); diff(u, y)]
% f = - diff(F(1), x) - diff(F(2), y)
% simplify(expand(f))
% g0y = subs(-Ku(1), x, 0)
% g1y = subs(Ku(1), x, 1)
% gx0 = subs(-Ku(2), y, 0)
% gx1 = subs(Ku(2), y, 1)

    function u = u(x, y)
        u = sin(2*pi*x) * exp(-2*pi/sqrt(delta)*y);
    end

    function du = du(~, ~)
        du = [2*pi*cos(2*pi*x) * exp(-2*pi/sqrt(delta)*y); ...
             -2*pi*sin(2*pi*x) * exp(-2*pi/sqrt(delta)*y) / sqrt(delta)];
    end

    function a = a(~, ~, ~)
        a = [1, 0; 0, delta];
    end

    function f = f(~, ~, ~)
        f = 0;
    end

PDE = struct('a', @a, 'u', @u, 'f', @f, 'du', @du);
end