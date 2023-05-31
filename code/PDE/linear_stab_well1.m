function PDE = linear_stab_well1()

inner_boundary = 12;
outer_boundary = 10;
theta = 3 / 8 * pi;
Kx = 1000;
G = [cos(theta), sin(theta); -sin(theta), cos(theta)];
kappa = G' * [Kx, 0; 0, 1] * G;

    % there is not exact solution
    function u = u(x, y)
        if x < 0.8 && x > 0.2 && ...
           y < 0.8 && y > 0.2
            u = inner_boundary;
        else
            u = outer_boundary;
        end
    end

    function du = du(~, ~)
        du = [0; 0];
    end

    function a = a(~, ~, ~)
        a = kappa;
    end

    function f = f(~, ~, ~)
        f = 0;
    end

PDE = struct('a', @a, 'u', @u, 'f', @f, 'du', @du, ...
    'bdtype', 'D');
end