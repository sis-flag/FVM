% function PDE = linear_stab5(a1, b1, a2, b2)
% 
% if nargin < 1
%     a1 = 1e2; b1 = 1e1; a2 = 1e-2; b2 = 1e-3;
% end
% 
% K1 = [a1, 0; 0, b1];
% K2 = [a2, 0; 0, b2];
% 
%     function u = u(x, ~)
%         % there is no analytical solution
%         u = 1 - x;
%     end
% 
%     function a = a(x, y, ~)
%         if x < 0.5
%             inty = floor(10 * (y + 0.15));
%             if mod(inty, 2) == 0
%                 a = K1;
%             else
%                 a = K2;
%             end
%         else
%             inty = floor(10 * y);
%             if mod(inty, 2) == 0
%                 a = K1;
%             else
%                 a = K2;
%             end
%         end
%     end
% 
%     function f = f(~, ~, ~)
%         f = 0;
%     end
% 
% 
% PDE = struct('a', @a, 'u', @u, 'f', @f);
% end