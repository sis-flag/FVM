function PDE = linear_stab10(delta)

if nargin < 1
    delta = 10;
end

ax = delta.^[1, -1, -2, 2];
ay = delta.^[-2, 2, 1, -1];
alpha = delta.^[-1 ,1, 2, -2];


    function u = u(x, y)
        if x < 0.5 && y < 0.5
            k = 1;
        elseif  x >= 0.5 && y < 0.5
            k = 2;
        elseif  x >= 0.5 && y >= 0.5
            k = 3;
        elseif  x < 0.5 && y >= 0.5
            k = 4;
        end
        u = alpha(k) * sin(2*pi*x) * sin(2*pi*y);
    end

    function du = du(x, y)
        if x < 0.5 && y < 0.5
            k = 1;
        elseif  x >= 0.5 && y < 0.5
            k = 2;
        elseif  x >= 0.5 && y >= 0.5
            k = 3;
        elseif  x < 0.5 && y >= 0.5
            k = 4;
        end
        du = alpha(k)*pi * [cos(pi*x)*sin(pi*y); cos(pi*y)*sin(pi*x)];
    end

    function a = a(x, y, ~)
        if x < 0.5 && y < 0.5
            k = 1;
        elseif  x >= 0.5 && y < 0.5
            k = 2;
        elseif  x >= 0.5 && y >= 0.5
            k = 3;
        elseif  x < 0.5 && y >= 0.5
            k = 4;
        end
        a = [ax(k), 0; 0, ay(k)];
    end

    function f = f(x, y, ~)
        if x < 0.5 && y < 0.5
            k = 1;
        elseif  x >= 0.5 && y < 0.5
            k = 2;
        elseif  x >= 0.5 && y >= 0.5
            k = 3;
        elseif  x < 0.5 && y >= 0.5
            k = 4;
        end
        f = 4*pi*pi * (ax(k)+ay(k)) * ...
            alpha(k) * sin(2*pi*x) * sin(2*pi*y);
    end


PDE = struct('a', @a, 'u', @u, 'f', @f, 'du', @du);
end