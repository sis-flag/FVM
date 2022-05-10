function PDE = problem1_2()

    function u = u(x, y)
        u = sin((x-1)*(y-1)) - (x-1)^3 * (y-1)^2;
    end

    function a = a(~, ~)
        a = [1.5, 0.5; 0.5, 1.5];
    end

    function f = f(x, y)
        f = (1.5*x^2 + x*y + 1.5*y^2 - 4*x - 4*y + 4) * sin((x-1)*(y-1)) ...
            + 30*x + 24*y - 30*x*y + 9*x*y^2 + 6*x^2*y - 15*x^2 + 3*x^3 - 9*y^2 - 18 ...
            - cos((x-1)*(y-1));
    end

PDE = struct('bdtype', 'D', 'a', @a, 'u', @u, 'f', @f);
end