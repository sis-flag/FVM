function PDE = linear_stab0()

    function u = u(x, y)
        u = x + y + 1;
    end

    function du = du(~, ~)
        du = [1; 1];
    end

    function a = a(~, ~, ~)
        a = [1, 0; 0, 1];
    end

    function f = f(~, ~, ~)
        f = 0;
    end


PDE = struct('a', @a, 'u', @u, 'f', @f, 'du', @du);
end