function PDE = linear_stab11()

aa = [2; 4; 4; 2];
bb = [6; 6; 8; 8];
cc = 0.5 * (aa + bb);
alpha = [1000; 500; 50; 100];
beta = [3; 2; 1; 2];
gamma = [1; 1; 2; 2];


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
        u = aa(k)*x + bb(k)*y - cc(k);
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
        du = [aa(k); bb(k)];
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
        a = [alpha(k), gamma(k); gamma(k), beta(k)];
    end

    function f = f(~, ~, ~)
        f = 0;
    end

PDE = struct('a', @a, 'u', @u, 'f', @f, 'du', @du, ...
    'bdtype', 'D');
end