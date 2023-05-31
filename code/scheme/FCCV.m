% function [Amat, F] = FCCV(Mesh, PDE, gamma)
% 
% if nargin < 3
%     gamma = 1;
% end
% 
% all_a = cell(Mesh.nU, 1);
% for U = 1:Mesh.nU
%     xc = Mesh.xc(U); yc = Mesh.yc(U);
%     all_a{U} = PDE.a(xc, yc, 0);
% end
% 
% F = zeros(Mesh.nP, 1);
% for P = 1:Mesh.nP
%     xp = Mesh.xp(P); yp = Mesh.yp(P);
%     
%     if Mesh.isbdp(P)
%         F(P) = PDE.u(xp, yp);
%         
%     else
%         nEs = Mesh.P2E{P};
%         nUs = Mesh.P2U{P};
%         
%         for U = nUs
%             tnE = intersect(nEs, Mesh.U2E{U});
%             
%             xe1 = Mesh.xe(tnE(1)); ye1 = Mesh.ye(tnE(1));
%             xe2 = Mesh.xe(tnE(2)); ye2 = Mesh.ye(tnE(2));
%             
%             xc = Mesh.xc(U); yc = Mesh.yc(U);
% 
%             xt1 = (xe1 + xp + xc) / 3;
%             yt1 = (ye1 + yp + yc) / 3;
%             Tarea1 = det([xe1-xp, xe1-xc; ye1-yp, ye1-yc]);
%             F1 = 0.5 * abs(Tarea1) * PDE.f(xt1, yt1, 0);
% 
%             xt2 = (xe2 + xp + xc) / 3;
%             yt2 = (ye2 + yp + yc) / 3;
%             Tarea2 = det([xe2-xp, xe2-xc; ye2-yp, ye2-yc]);
%             F2 = 0.5 * abs(Tarea2) * PDE.f(xt2, yt2, 0);
% 
%             F(P) = F(P) + F1 + F2;
%         end
%     end
% end
% 
% % Amat = spalloc(Mesh.nP, Mesh.nP, 50*Mesh.nP);
% iA = zeros(50*Mesh.nP, 1);
% jA = zeros(50*Mesh.nP, 1);
% vA = zeros(50*Mesh.nP, 1);
% nnzA = 0;
% for U = 1:Mesh.nU
%     nE = Mesh.U2E{U};
%     nP = Mesh.U2P{U};
%     
%     Kinv = 1 / Mesh.area(U);
%     xP = Mesh.xp(nP); yP = Mesh.yp(nP);
%     xE = Mesh.xe(nE); yE = Mesh.ye(nE);
%     xU = Mesh.xc(U); yU = Mesh.yc(U);
%     
%     n = length(nP);
%     R = -eye(n) + diag(ones(n-1,1), 1);
%     R(n, 1) = 1;
%     
%     X = R * [xP, yP];
%     N = [xE - xU, yE - yU] * [0,1;-1,0];
%     C = eye(n) - X / (X'*X) * X';
%     A = Kinv * N * all_a{U} * N' + gamma * (C * C');
%     RAR = R' * A * R;
%     
%     for r = 1:length(nP)
%         if ~Mesh.isbdp(nP(r))
%             for c = 1:length(nP)
%                 nnzA = nnzA + 1;
%                 iA(nnzA) = nP(r); jA(nnzA) = nP(c);
%                 vA(nnzA) = RAR(r, c);
%             end
%         end
%     end
% end
% 
% for P = 1:Mesh.nP
%     if Mesh.isbdp(P)
%         nnzA = nnzA + 1;
%         iA(nnzA) = P; jA(nnzA) = P;
%         vA(nnzA) = 1;
%     end
% end
% 
% Amat = sparse(iA(1:nnzA), jA(1:nnzA), vA(1:nnzA), ...
%     Mesh.nP, Mesh.nP);
% 
% end