function PDE = linear_stab0()

    function u = u(x, y)
        u = 3*x + 2*y + 1;
    end

    function du = du(~, ~)
        du = [3; 2];
    end

    function a = a(~, ~, ~)
        a = [1.5, 0.5; 0.5, 1.5];
    end

    function f = f(~, ~, ~)
        f = 0;
    end

PDE = struct('a', @a, 'u', @u, 'f', @f, 'du', @du, ...
    'bdtype', 'D');
end