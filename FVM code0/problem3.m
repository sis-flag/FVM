function PDE = problem3(delta, theta)

if nargin < 1
    delta = 1e-3;
end
if nargin < 2
    theta = 2/9*pi;
end

R = [cos(theta), sin(theta); -sin(theta), cos(theta)];

    function a = a(~, ~)
        a = R' * [1, 0; 0, delta] * R;
    end

    function u = u(x, y)
        % there is no analytical solution
        if x == 0
            if y <= 0.2 % y in [0, 0.2]
                u = 1;
            elseif y <= 0.3 % y in [0.2, 0.3]
                u = 2 - 5 * y;
            else % y in [0.3, 1]
                u = 0.5;
            end
        elseif x == 1
            if y <= 0.7
                u = 0.5;
            elseif y <= 0.8
                u = 4 - 5 * y;
            else
                u = 0;
            end
        elseif y == 0
            if x <= 0.2
                u = 1;
            elseif x <= 0.3
                u = 2 - 5 * x;
            else
                u = 0.5;
            end
        elseif y == 1
            if x <= 0.7
                u = 0.5;
            elseif x <= 0.8
                u = 4 - 5 * x;
            else
                u = 0;
            end
        end
    end

    function f = f(~, ~)
        f = 0;
    end


PDE = struct('bdtype', 'D', 'a', @a, 'u', @u, 'f', @f);
end