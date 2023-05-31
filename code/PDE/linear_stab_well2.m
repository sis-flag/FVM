function PDE = linear_stab_well2()

left_boundary = 12;
right_boundary = 10;
theta = 3 / 8 * pi;
Kx = 1000;
G = [cos(theta), sin(theta); -sin(theta), cos(theta)];
kappa = G' * [Kx, 0; 0, 1] * G;

    % there is no exact solution
    function u = u(x, ~)
        if x < 0.5
            u = left_boundary;
        else
            u = right_boundary;
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
    'bdtype', 'N');
end