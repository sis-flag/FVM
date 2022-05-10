function phi = mv2dglob(v,x)
%
% Evaluate MV basis functions in a convex polygon
% using a global formula that is valid also on the boundary.
%
% Inputs:
% v : [x1 y1; x2 y2; ...; xn yn], vertices of the polygon, ccw
% x : [x(1) x(2)], point inside polygon
% Outputs:
% phi : basis functions = [phi_1; ...; phi_n]46 M. S. Floater
n = size(v,1);
w = zeros(n,1);
phi = zeros(n,1);
r = zeros(n,1);
d = zeros(n,2);
for i = 1:n
    d(i,:) = v(i,:) - x;
    r(i) = norm(d(i,:));
end
mod1 = inline('mod(i-1,n)+1','i','n');
for i = 1:n
    im1 = mod1(i-1,n);
    ip1 = mod1(i+1,n);
    prod = 1;
    for j=1:n
        if j ~= im1 && j ~= i
            jp1 = mod1(j+1,n);
            prod = prod * (r(j) * r(jp1) + dot(d(j,:),d(jp1,:)));
        end
    end
    prod = prod * (r(im1) * r(ip1) - dot(d(im1,:),d(ip1,:)));
    w(i) = sqrt(prod);
end
wsum = sum(w);
phi = w/wsum;