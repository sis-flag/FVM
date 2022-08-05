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

    function E0 = find_edge(P1, P2)
        % find the edge connecting P1 and P2
        % if there is no such edge, return -1
        for E1 = P2E{P1}
            nP1 = E2P{E1};
            if (nP1(1) == P1 && nP1(2) == P2) || ...
                    (nP1(1) == P2 && nP1(2) == P1)
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
            P1 = nPs(end);
            P2 = nPs(1);
        else
            P1 = nPs(k-1);
            P2 = nPs(k);
        end
        
        E = find_edge(P1, P2);
        
        if E == -1
            % the edge is not recorded yet, generate a new edge
            nE = nE + 1;
            E = nE;
            E2P{E}= [P1, P2];
            P2E{P1}(end+1) = E;
            P2E{P2}(end+1) = E;
        end
        
        if length(E2U) < E
            E2U{E} = [];
        end
        E2U{E}(end+1) = U;
        U2E{U}(end+1) = E;
    end
end

xp = zeros(nP, 1); % coordinats of points
yp = zeros(nP, 1); % coordinats of points
xe = zeros(nE, 1); % coordinats of the middle points of edges
ye = zeros(nE, 1); % coordinats of the middle points of edges
xc = zeros(nU, 1); % coordinats of the center points of units
yc = zeros(nU, 1); % coordinats of the center points of units
nx = zeros(nE, 1); % normal vectors of edges (point out from E2U{1})
ny = zeros(nU, 1); % normal vectors of edges (point out from E2U{1})
len = zeros(nE, 1); % length of edges
area = zeros(nU, 1); % area of units
isbdp = zeros(nP, 1, 'logical'); % is (or not) boundary point
isbdu = zeros(nU, 1, 'logical'); % is (or not) boundary unit

for P = 1:nP
    xp(P) = coord(1, P);
    yp(P) = coord(2, P);
end

for U = 1:nU
    xc(U) = mean(coord(1, U2P{U}));
    yc(U) = mean(coord(2, U2P{U}));
    
    % area is calculated by the formula in
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
    
    area(U) = 0.5 * abs(tarea);
end

for E = 1:nE
    nPs = E2P{E};
    P1 = nPs(1); P2 = nPs(2);
    
    xe(E) = (coord(1, P2) + coord(1, P1)) / 2;
    ye(E) = (coord(2, P2) + coord(2, P1)) / 2;
    
    vec = coord(:, P2) - coord(:, P1);
    len(E) = norm(vec);
    nx(E) = vec(2) / len(E);
    ny(E) = -vec(1) / len(E);
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
Mesh.yp = yp;
Mesh.xe = xe;
Mesh.ye = ye;
Mesh.xc = xc;
Mesh.yc = yc;
Mesh.nx = nx;
Mesh.ny = ny;
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
