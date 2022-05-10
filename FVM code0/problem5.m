function PDE = problem5(delta)

if nargin < 1
    delta = 1e-3;
end

    function a = a(x, y)
        rt = x*x+y*y ;
        a = [delta*x*x+y*y,  -(1-delta)*x*y ; -(1-delta)*x*y, x*x+delta*y*y] / rt;
    end

    function u = u(x, y)
        u = sin(pi*x) * sin(pi*y);
    end

    function f = f(x, y)
        
        f = ...
            + x^4*pi^2*sin(pi*x)*sin(pi*y) ...
            + y^4*pi^2*sin(pi*x)*sin(pi*y) ...
            + x^3*pi*cos(pi*x)*sin(pi*y) ...
            + y^3*pi*cos(pi*y)*sin(pi*x) ...
            - delta*x^3*pi*cos(pi*x)*sin(pi*y) ...
            - delta*y^3*pi*cos(pi*y)*sin(pi*x) ...
            + x*y^2*pi*cos(pi*x)*sin(pi*y) ...
            + x^2*y*pi*cos(pi*y)*sin(pi*x) ...
            + 2*x*y^3*pi^2*cos(pi*x)*cos(pi*y) ...
            + 2*x^3*y*pi^2*cos(pi*x)*cos(pi*y) ...
            + delta*x^4*pi^2*sin(pi*x)*sin(pi*y) ...
            + delta*y^4*pi^2*sin(pi*x)*sin(pi*y) ...
            + 2*x^2*y^2*pi^2*sin(pi*x)*sin(pi*y) ...
            + 2*delta*x^2*y^2*pi^2*sin(pi*x)*sin(pi*y) ...
            - delta*x*y^2*pi*cos(pi*x)*sin(pi*y) ...
            - delta*x^2*y*pi*cos(pi*y)*sin(pi*x) ...
            - 2*delta*x*y^3*pi^2*cos(pi*x)*cos(pi*y) ...
            - 2*delta*x^3*y*pi^2*cos(pi*x)*cos(pi*y);
        f = f / (x*x + y*y)^2;
    end


PDE = struct('bdtype', 'D', 'a', @a, 'u', @u, 'f', @f);
end