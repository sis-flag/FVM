function Mesh = arrange_polygonal(coord, U2P)
% generate a mesh struct for polygonal mesh
% input:
%     coord (2 x n array): coordinates of each points.
%     U2P (size n cell): index of each neighbor points of each unit,
%     U2P must be ordered clockwise
% output: a mesh struct

nP = size(coord, 2); % the number of points
nU = length(U2P); % the number of units
nE = 0; % the number of edges

P2U = cell(nP, 1); % neighbor units of each point
P2E = cell(nP, 1); % neighbor edges of each point
E2P = cell(nE, 1); % neighbor points of each edge
E2U = cell(nE, 1); % neighbor units of each edge
U2E = cell(nU, 1); % neighbor edges of each unit

for U = 1:nU
    for P = U2P{U}
        P2U{P}(end+1) = U;
    end
end

    function E0 = find_edge(P0, P1)
        % find the edge connecting P0 and P1
        % if there is no such edge, return -1
        for E1 = P2E{P0}
            nP1 = E2P{E1};
            if (nP1(1) == P0 && nP1(2) == P1) || ...
                    (nP1(1) == P1 && nP1(2) == P0)
                E0 = E1;
                return
            end
        end
        E0 = -1;
    end

for U = 1:nU
    nPs = U2P{U};
    
    for k = 1:length(nPs)
        % all edges surrounding U
        if k == 1
            P0 = nPs(end);
            P1 = nPs(1);
        else
            P0 = nPs(k-1);
            P1 = nPs(k);
        end
        
        E = find_edge(P0, P1);
        
        if E == -1
            % the edge is not recorded yet, generate a new edge
            nE = nE + 1;
            E = nE;
            E2P{E}= [P0, P1];
            P2E{P0}(end+1) = E;
            P2E{P1}(end+1) = E;
        end
        
        if length(E2U) < E
            E2U{E} = [];
        end
        E2U{E}(end+1) = U;
        U2E{U}(end+1) = E;
    end
end

xp = cell(nP, 1); % coordinats of points
xe = cell(nE, 1); % coordinats of the middle points of edges
xc = cell(nU, 1); % coordinats of the center points of units
nv = cell(nE, 1); % normal vectors of edges (point out from E2U{1})
len = cell(nE, 1); % length of edges
area = cell(nU, 1); % area of units
isbdp = zeros(nP, 1, 'logical'); % is (or not) boundary point
isbdu = zeros(nU, 1, 'logical'); % is (or not) boundary unit

for P = 1:nP
    xp{P} = coord(:, P);
end

for U = 1:nU
    xc{U} = mean(coord(:,U2P{U}), 2);
    
    % area is calculated by the formula here
    % https://zhuanlan.zhihu.com/p/110025234
    nPs = U2P{U};
    tarea = 0;
    
    for k = 1:length(nPs)
        if k == 1
            xp0 = coord(:, nPs(end));
            xp1 = coord(:, nPs(1));
        else
            xp0 = coord(:, nPs(k-1));
            xp1 = coord(:, nPs(k));
        end
        tarea = tarea + det([xp0, xp1]);
    end
    
    area{U} = 0.5 * abs(tarea);
end

for E = 1:nE
    nPs = E2P{E};
    xe{E} = (coord(:,nPs(2)) + coord(:,nPs(1))) / 2;
    vec = coord(:,nPs(2)) - coord(:,nPs(1));
    len{E} = norm(vec);
    nv{E} = [0, 1; -1, 0] * vec / len{E};
end

for E = 1:nE
    if length(E2U{E}) == 1
        isbdu(E2U{E}) = 1;
        for P = E2P{E}
            isbdp(E2P{E}) = 1;
        end
    end
end

Mesh = struct();
Mesh.nP = nP;
Mesh.nU = nU;
Mesh.nE = nE;
Mesh.xp = xp;
Mesh.xe = xe;
Mesh.xc = xc;
Mesh.nv = nv;
Mesh.len = len;
Mesh.area = area;
Mesh.isbdp = isbdp;
Mesh.isbdu = isbdu;
Mesh.P2U = P2U;
Mesh.P2E = P2E;
Mesh.E2P = E2P;
Mesh.E2U = E2U;
Mesh.U2E = U2E;
Mesh.U2P = U2P;

end
