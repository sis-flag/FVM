function PDE = problem6(delta, a1, b1, a2, b2)

if nargin < 1
    delta = 0.2;
    a1 = 1; b1 = 0.1; a2 = 100; b2 = 10;
end

cost = 1 / sqrt(1 + delta*delta);
sint = delta * cost;
R = [cost, sint; -sint, cost];

K1 = R' * [a1, 0; 0, b1] * R;
K2 = R' * [a2, 0; 0, b2] * R;

phi1 = @(x,y) y - delta*(x-0.5) - 0.475;
phi2 = @(x,y) phi1(x,y) - 0.05;

    function a = a(x, y)
        if phi1(x,y) < 0 || phi2(x,y) > 0
            a = K1;
        else
            a = K2;
        end
    end

    function u = u(x, y)
        u = - x - y * delta;
    end

    function f = f(~, ~)
        f = 0;
    end

PDE = struct('bdtype', 'D', 'a', @a, 'u', @u, 'f', @f);
end