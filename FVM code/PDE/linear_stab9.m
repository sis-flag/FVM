function PDE = linear_stab9(delta, fn)

if nargin < 1
    delta = 4;
end
if nargin < 2
    fn = 2;
end

    function u = u(x, y)
        if x < 0.5
            u = (4+4*y-fn*y*y) + (2*x-1)*(3+2*y);
        else
            u = (4+4*y-fn*y*y) + (2*x-1)*(3+2*y)/delta;
        end
    end

    function du = du(x, y)
        if x < 0.5
            du = [4*y+6; 4*x-2*fn*y+2];
        else
            du = [0; 4-2*fn*y] + [4*y+6; 4*x-2]/delta;
        end
    end

    function a = a(x, ~, ~)
        if x < 0.5
            a = [1, 0; 0, 1];
        else
            a = delta * [1, 0; 0, 1];
        end
    end

    function f = f(x, ~, ~)
        if x < 0.5
            f = 2* fn;
        else
            f = 2* delta * fn;
        end
    end


PDE = struct('a', @a, 'u', @u, 'f', @f, 'du', @du);
end