function Mesh = get_circle_mesh(Nr, Nt)

r = linspace(0, 1, Nr+1);
t = linspace(0, 2*pi, Nt+1);
[r, t] = meshgrid(r(2:end), t(1:end-1));

x = r.*cos(t); y = r.*sin(t);

coord = [0, x(:)'; 0, y(:)'];

U2P = cell(Nr*Nt, 1);
for nr = 1:Nr
    for nt = 1:Nt
        k = (nr-1)*Nt + nt;
        if nr == 1
            if nt < Nt
                U2P{k} = [1, nt+1, nt+2];
            else
                U2P{k} = [1, Nt+1, 2];
            end
        else
            if nt < Nt
                U2P{k} = [k-Nt+2, k-Nt+1, k+1, k+2];
            else
                U2P{k} = [k-Nt-Nt+2, k-Nt+1, k+1, k+2-Nt];
            end
        end
    end
end


Mesh = arrange_polygonal(coord, U2P);

end