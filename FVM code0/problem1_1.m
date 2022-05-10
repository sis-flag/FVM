function PDE = problem1_1()

% syms x y
% u = 16 * x * (1-x) * y * (1-y)
% a = [1.5, 0.5; 0.5, 1.5]
% F = a * [diff(u, x); diff(u, y)]
% f = - diff(F(1), x) - diff(F(2), y)
% simplify(expand(f))

    function u = u(x, y)
        u = 16 * x * (1-x) * y * (1-y);
    end

    function a = a(~, ~)
        a = [1.5, 0.5; 0.5, 1.5];
    end

    function f = f(x, y)
        f = - 48*x^2 - 64*x*y + 80*x - 48*y^2 + 80*y - 16;
    end

PDE = struct('bdtype', 'D', 'a', @a, 'u', @u, 'f', @f);
end
