function PDE = linear_stab3(delta)

if nargin < 1
    delta = 1e5;
end

% syms x y delta
% u = sin(2*pi*x) * exp(-2*pi/sqrt(delta)*y)
% K = [1, 0; 0, delta]
% Ku = K * [diff(u, x); diff(u, y)]
% f = - diff(Ku(1), x) - diff(Ku(2), y)
% simplify(expand(f))
% g0y = subs(-Ku(1), x, 0)
% g1y = subs(Ku(1), x, 1)
% gx0 = subs(-Ku(2), y, 0)
% gx1 = subs(Ku(2), y, 1)

    function u = u(x, y)
        u = sin(2*pi*x) * exp(-2*pi/sqrt(delta)*y);
    end

    function a = a(~, ~, ~)
        a = [1, 0; 0, delta];
    end

    function f = f(~, ~, ~)
        f = 0;
    end

%     function g = g(x, y)
%         if x == 0
%             g = -2*pi*exp(-(2*pi*y)/delta^(1/2));
%         elseif x == 1
%             g = 2*pi*exp(-(2*pi*y)/delta^(1/2));
%         elseif y == 0
%             g = 2*delta^(1/2)*pi*sin(2*pi*x);
%         elseif y == 1
%             g = -2*delta^(1/2)*pi*exp(-(2*pi)/delta^(1/2))*sin(2*pi*x);
%         end
%     end

PDE = struct('a', @a, 'u', @u, 'f', @f);
end