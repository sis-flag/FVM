% function weight = interp_my_weight(Mesh, PDE, uc)
% % get interplation weights on mesh points
% 
% if nargin < 3
%     uc = zeros(Mesh.nU, 1);
% end
% 
% all_a = cell(Mesh.nU, 1);
% for U = 1:Mesh.nU
%     xu = Mesh.xc(U); yu = Mesh.yc(U);
%     all_a{U} = PDE.a(xu, yu, uc(U));
% end
% 
% weight = cell(Mesh.nP, 1);
% for P = 1:Mesh.nP
%     if Mesh.isbdp(P)
%         weight{P} = [];
%     else
%         nU = Mesh.P2U{P};
%         xp = Mesh.xp(P); yp = Mesh.yp(P);
%         
%         dist_inv = zeros(length(nU), 1);
%         for k = 1:length(nU)
%             xU = Mesh.xc(nU(k));
%             dist_inv(k) = 1 / norm(xU - xP);
%         end
%         w0 = dist_inv / sum(dist_inv);
%         
%         A = zeros(3, length(nU));
%         A(1, :) = 1;
%         for k = 1:length(nU)
%             xu = Mesh.xc(nU(k));
%             yu = Mesh.yc(nU(k));
%             A(2, k) = xu - xp;
%             A(3, k) = yu - yp;
%         end
%         
%         b = [1; 0; 0];
%         
%         lambda = (A * A') \ (A * w0 - b);
%         w = w0 - A' * lambda;
%         
%         if any(w < 0)
%             weight{P} = w0;
%         else
%             weight{P} = w;
%         end
%     end
% end
% 
% end