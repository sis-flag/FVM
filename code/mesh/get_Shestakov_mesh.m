function Mesh = get_Shestakov_mesh(N, p)
% generate a Shestakov mesh with size 2^N on domain [0,1] x [0,1]

% This routine generates an nm by nm mesh, commonly known as a
% Shestakov mesh, with 0 < r < rm and 0 < z < zm. The randomness
% is controlled by parameter "a", where 0 <= a <= 0.5. A value of 
% p = 0.5 gives a rectangular mesh. Boomerang zones are allowed 
% but bowties are not allowed.
%
% Matlab version by Jean C. Ragusa, Texas A&M University, April 1, 2013
%
% based on the original Fortran program by Alek Shestakov; that program 
% is found in: Nuclear Science & Engineering, Vol 105, pp.88-104 (1990),
% "Test Problems in Radiative Transfer Calculations", by A. Shestakov, 
% D. Kershaw, and G. Zimmerman

% input:
%     N (integer): number of subdivisions of the original rectangle
%     p (nummber): perturbation of the mesh (0 < p < 0.5)
% output:
%     mesh with size (2^N+1, 2^N+1)

% allocate memory
Nm = 2^N + 1;
ind = Nm - 1;
x=zeros(Nm,Nm);
y=zeros(Nm,Nm);

% initialize 4 corners
for i = 0:1
    k = 1 + i*ind;
    for j = 0:1
        l = 1 + j*ind;
        x(k,l) = (1-i);
        y(k,l) = j;
    end
end
% fill in the rest of the points
for nl = 0:N-1
    nn = 2^nl;
    inc = ind/nn;
    inch = inc/2;
    for k = 1:nn
        k1 = 1+(k-1)*inc;
        k2 = k1+inc;
        k3 = k1+inch;
        for l = 1:nn
            l1 = 1+(l-1)*inc;
            l2 = l1+inc;
            l3 = l1+inch;
            if (l==1)
                ar = p+rand(1,1)*(1-2*p);
                x(k3,l1) = ar*x(k1,l1)+(1-ar)*x(k2,l1);
                y(k3,l1) = ar*y(k1,l1)+(1-ar)*y(k2,l1);
            end
            ar = p+rand(1,1)*(1-2*p);
            x(k2,l3) = ar*x(k2,l1)+(1-ar)*x(k2,l2);
            y(k2,l3) = ar*y(k2,l1)+(1-ar)*y(k2,l2);
            ar = p+rand(1,1)*(1-2*p);
            x(k3,l2) = ar*x(k1,l2)+(1-ar)*x(k2,l2);
            y(k3,l2) = ar*y(k1,l2)+(1-ar)*y(k2,l2);
            if (k==1)
                ar = p+rand(1,1)*(1-2*p);
                x(k1,l3) = ar*x(k1,l1)+(1-ar)*x(k1,l2);
                y(k1,l3) = ar*y(k1,l1)+(1-ar)*y(k1,l2);
            end
            ar = p+rand(1,1)*(1-2*p);
            br = p+rand(1,1)*(1-2*p);
            r1 = x(k1,l1);
            r2 = x(k2,l1);
            r3 = x(k2,l2);
            r4 = x(k1,l2);
            z1 = y(k1,l1);
            z2 = y(k2,l1);
            z3 = y(k2,l2);
            z4 = y(k1,l2);
            % check for boomerang zones
            det2 = (r2-r1)*(z3-z2)-(r3-r2)*(z2-z1);
            det3 = (r3-r2)*(z4-z3)-(r4-r3)*(z3-z2);
            det4 = (r4-r3)*(z1-z4)-(r1-r4)*(z4-z3);
            det1 = (r1-r4)*(z2-z1)-(r2-r1)*(z1-z4);
            if (det2>0)
                d = (r4-r3)*(z2-z1)-(r2-r1)*(z4-z3);
                r3p = ((r2-r1)*(r4*z3-r3*z4)-(r4-r3)*(r2*z1-r1*z2))/d;
                z3p = ((z2-z1)*(r4*z3-r3*z4)-(z4-z3)*(r2*z1-r1*z2))/d;
                d = (r4-r1)*(z2-z3)-(r2-r3)*(z4-z1);
                r1p = ((r2-r3)*(r4*z1-r1*z4)-(r4-r1)*(r2*z3-r3*z2))/d;
                z1p = ((z2-z3)*(r4*z1-r1*z4)-(z4-z1)*(r2*z3-r3*z2))/d;
                r3 = r3p;
                z3 = z3p;
                r1 = r1p;
                z1 = z1p;
            elseif (det3>0)
                d = (r1-r4)*(z3-z2)-(r3-r2)*(z1-z4);
                r4p = ((r3-r2)*(r1*z4-r4*z1)-(r1-r4)*(r3*z2-r2*z3))/d;
                z4p = ((z3-z2)*(r1*z4-r4*z1)-(z1-z4)*(r3*z2-r2*z3))/d;
                d = (r1-r2)*(z3-z4)-(r3-r4)*(z1-z2);
                r2p = ((r3-r4)*(r1*z2-r2*z1)-(r1-r2)*(r3*z4-r4*z3))/d;
                z2p = ((z3-z4)*(r1*z2-r2*z1)-(z1-z2)*(r3*z4-r4*z3))/d;
                r4 = r4p;
                z4 = z4p;
                r2 = r2p;
                z2 = z2p;
            elseif (det4>0)
                d = (r2-r1)*(z4-z3)-(r4-r3)*(z2-z1);
                r1p = ((r4-r3)*(r2*z1-r1*z2)-(r2-r1)*(r4*z3-r3*z4))/d;
                z1p = ((z4-z3)*(r2*z1-r1*z2)-(z2-z1)*(r4*z3-r3*z4))/d;
                d = (r2-r3)*(z4-z1)-(r4-r1)*(z2-z3);
                r3p = ((r4-r1)*(r2*z3-r3*z2)-(r2-r3)*(r4*z1-r1*z4))/d;
                z3p = ((z4-z1)*(r2*z3-r3*z2)-(z2-z3)*(r4*z1-r1*z4))/d;
                r1 = r1p;
                z1 = z1p;
                r3 = r3p;
                z3 = z3p;
            elseif (det1>0)
                d = (r3-r2)*(z1-z4)-(r1-r4)*(z3-z2);
                r2p = ((r1-r4)*(r3*z2-r2*z3)-(r3-r2)*(r1*z4-r4*z1))/d;
                z2p = ((z1-z4)*(r3*z2-r2*z3)-(z3-z2)*(r1*z4-r4*z1))/d;
                d = (r3-r4)*(z1-z2)-(r1-r2)*(z3-z4);
                r4p = ((r1-r2)*(r3*z4-r4*z3)-(r3-r4)*(r1*z2-r2*z1))/d;
                z4p = ((z1-z2)*(r3*z4-r4*z3)-(z3-z4)*(r1*z2-r2*z1))/d;
                r2 = r2p;
                z2 = z2p;
                r4 = r4p;
                z4 = z4p;
            end
            x(k3,l3) = ar*br*r1 + ar*(1-br)*r4 + (1-ar)*br*r2 + (1-ar)*(1-br)*r3;
            y(k3,l3) = ar*br*z1 + ar*(1-br)*z4 + (1-ar)*br*z2 + (1-ar)*(1-br)*z3;
        end
    end
end

Mesh = arrange_quadrilateral(x(end:-1:1, :), y(end:-1:1, :));
end