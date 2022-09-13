% function PDE = linear_stab8(delta, alpha)
% 
% if nargin < 2
%      alpha = 1e-2;
% end
% if nargin < 1
%     delta = 0.2;
% end
% 
% phi1 = @(x,y) y - delta*(x-0.5) - 0.475;
% phi2 = @(x,y) phi1(x,y) - 0.05;
% 
%     function u = u(x, y)
%         if  phi1(x,y) < 0
%             u = - phi1(x,y);
%         elseif phi2(x,y) < 0
%             u = - phi1(x,y)/0.01;
%         else
%             u = - phi2(x,y) - 5;
%         end
%     end
% 
%     function a = a(x, y, ~)
%         if  phi1(x,y) < 0
%             a = [1, 0; 0, 1];
%         elseif phi2(x,y) < 0
%             a = [alpha, 0; 0, alpha];
%         else
%             a = [1, 0; 0, 1];
%         end
%     end
% 
%     function f = f(~, ~, ~)
%         f = 0;
%     end
% 
% PDE = struct('a', @a, 'u', @u, 'f', @f);
% end