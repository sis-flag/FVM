clear;

%% equation 1-1
e = 0.1;
b = 0.1;
k = 1;
ue = @(x) sin(pi*x) + b;
a = @(x,u) k * u + e;
au = @(x,u) k;
f = @(x,u) pi*pi*((e+b*k)*sin(pi*x) + 2*k*sin(pi*x)^2 - k);
fu = @(x,u) 0;

%% equation 1-2
% e = 0.1;
% b = 0.1;
% k = 1;
% ue = @(x) sin(pi*x) + b;
% a = @(x,u) k * u + e;
% au = @(x,u) k;
% f = @(x,u) pi*pi*((2*k*u-k*b+e) * sin(pi*x) - k);
% fu = @(x,u) 2*pi*pi*k*sin(pi*x);

%% equation 2-1
% e = 0.1;
% b = 0.1;
% ue = @(x) sin(pi*x*x) + b;
% a = @(x,u) exp(u) - 1 + e;
% au = @(x,u) exp(u);
% f = @(x,u) 4*x^2*pi^2*sin(pi*x^2)*(e+exp(b+sin(pi*x^2))-1) ...
%     - 4*x^2*pi^2*exp(b + sin(pi*x^2))*cos(pi*x^2)^2 ...
%     - 2*pi*cos(pi*x^2)*(e + exp(b + sin(pi*x^2)) - 1);
% fu = @(x,u) 0;

%% equation 2-2
% e = 0.1;
% b = 1;
% ue = @(x) sin(pi*x*x) + b;
% a = @(x,u) exp(u) - 1 + e;
% au = @(x,u) exp(u);
% f = @(x,u) 4*x^2*pi^2*exp(u)*((b - u)^2 - 1) ...
%     - 2*pi*cos(pi*x^2)*(e + exp(u) - 1) ...
%     + 4*x^2*pi^2*sin(pi*x^2)*(e + exp(u) - 1);
% fu = @(x,u) 4*x^2*pi^2*exp(u)*((b - u)^2 - 1) ...
%     - 2*pi*exp(u)*cos(pi*x^2) ...
%     - 4*x^2*pi^2*exp(u)*(2*b - 2*u) ...
%     + 4*x^2*pi^2*exp(u)*sin(pi*x^2);

%% mesh
N = 50;
h = 1 / N;
x = linspace(0, 1, N+1)';

%% nolinear functions
A = @(U) getA(a, x, [ue(0); U; ue(1)]);
F = @(U) getF(f, a, x, [ue(0); U; ue(1)]);
Phi = @(U) A(U) * U - F(U);
J = @(U) getJ(a, au, fu, x, [ue(0); U; ue(1)]);

%% Picard
U0 = ones(N-1, 1);

aU = zeros(N-1, 50+1, 'double');
aP = zeros(N-1, 50+1, 'double');
aU(:,1) = U0;
aP(:,1) = Phi(U0);
for s = 1:50
    
    U1 = A(U0) \ F(U0);
    
    aU(:, s+1) = U1;
    aP(:, s+1) = Phi(U1);
    if sqrt(h)*norm(Phi(U1)) < 1e-12 ||...
            sqrt(h)*norm(U0-U1) < 1e-12
        break
    end
    U0 = U1;
end

Nit_P = s;
Err_P = sqrt(h*sum((aU(:,1:s) - aU(:,s+1)).^2));
Phi_P = sqrt(h*sum(aP(:,1:s+1).^2));

disp('Picard:')
disp(s)

%% Stablized Picard
U0 = ones(N-1, 1);
delta = 0.01;

aU = zeros(N-1, 50+1, 'double');
aP = zeros(N-1, 50+1, 'double');
aU(:,1) = U0;
aP(:,1) = Phi(U0);
for s = 1:50
    
    U1 = (A(U0) + delta * eye(N-1)) \ (F(U0) + delta * U0);
    
    aU(:, s+1) = U1;
    aP(:, s+1) = Phi(U1);
    if sqrt(h)*norm(Phi(U1)) < 1e-12 ||...
            sqrt(h)*norm(U0-U1) < 1e-12
        break
    end
    U0 = U1;
end

Nit_SP = s;
Err_SP = sqrt(h*sum((aU(:,1:s) - aU(:,s+1)).^2));
Phi_SP = sqrt(h*sum(aP(:,1:s+1).^2));

disp('Stablized Picard:')
disp(s)

%% Newton
U0 = ones(N-1, 1);

aU = zeros(N-1, 50+1, 'double');
aP = zeros(N-1, 50+1, 'double');
aU(:,1) = U0;
aP(:,1) = Phi(U0);
for s = 1:50
    
    U1 = U0 - J(U0) \ Phi(U0);
    
    aU(:, s+1) = U1;
    aP(:, s+1) = Phi(U1);
    if sqrt(h)*norm(Phi(U1)) < 1e-12 ||...
            sqrt(h)*norm(U0-U1) < 1e-12
        break
    end
    U0 = U1;
end

Nit_N = s;
Err_N = sqrt(h*sum((aU(:,1:s) - aU(:,s+1)).^2));
Phi_N = sqrt(h*sum(aP(:,1:s+1).^2));

disp('Newton:')
disp(s)

%% Picard-Newton
% TODO

%% Broyden
U0 = ones(N-1, 1);
B0 = inv(A(U0));

aU = zeros(N-1, 50+1, 'double');
aP = zeros(N-1, 50+1, 'double');
aU(:,1) = U0;
aP(:,1) = Phi(U0);
for s = 1:50
    
    U1 = U0 - B0 * Phi(U0);
    dU = U1 - U0;
    dP = Phi(U1) - Phi(U0);
    B1 = B0 + (dU-B0*dP) * (dU'*B0) / (dU'*B0*dP);
    
    aU(:, s+1) = U1;
    aP(:, s+1) = Phi(U1);
    if sqrt(h)*norm(Phi(U1)) < 1e-12 ||...
            sqrt(h)*norm(U0-U1) < 1e-12
        break
    end
    U0 = U1;
    B0 = B1;
end

Nit_B = s;
Err_B = sqrt(h*sum((aU(:,1:s) - aU(:,s+1)).^2));
Phi_B = sqrt(h*sum(aP(:,1:s+1).^2));

disp('Broyden:')
disp(s)

%% JFNK
% TODO

%% Sparse Newton
U0 = ones(N-1, 1);
epsilon = 1e-9;
deepth = 4;

aU = zeros(N-1, 50+1, 'double');
aP = zeros(N-1, 50+1, 'double');
ErrJ_SN = zeros(1, 50+1, 'double');

aU(:,1) = U0;
aP(:,1) = Phi(U0);
J0 = A(U0);
ErrJ_SN(1) = norm(J0 - J(U0), 'fro');
for s = 1:50
    
    dU = zeros(N-1, deepth, 'double');
    for m = 1:deepth
        dU(:, m) = epsilon * rand(N-1, 1);
        dP(:, m) = Phi(U0 + dU(:,m)) - Phi(U0);
    end
    J0 = approxJ(dU, dP);
    
    U1 = U0 - J0 \ Phi(U0);
    
    aU(:, s+1) = U1;
    aP(:, s+1) = Phi(U1);
    ErrJ_SN(s+1) = norm(J0 - J(U0), 'fro');
    if sqrt(h)*norm(Phi(U1)) < 1e-12 ||...
            sqrt(h)*norm(U0-U1) < 1e-12
        break
    end
    U0 = U1;
end

Nit_SN = s;
Err_SN = sqrt(h*sum((aU(:,1:s) - aU(:,s+1)).^2));
Phi_SN = sqrt(h*sum(aP(:,1:s+1).^2));
ErrJ_SN = ErrJ_SN(1:s+1);

disp('Sparse Newton:')
disp(s)

%% Sparse Quasi-Newton
U0 = ones(N-1, 1);

aU = zeros(N-1, 50+1, 'double');
aP = zeros(N-1, 50+1, 'double');
ErrJ_SQN = zeros(1, 50+1, 'double');

aU(:,1) = U0;
aP(:,1) = Phi(U0);
J0 = A(U0);
ErrJ_SQN(1) = norm(J0 - J(U0), 'fro');
for s = 1:50
    
    if s < 5
        J0 = A(U0);
    else
        dU = diff(aU(:,s-4:s), 1, 2);
        dP = diff(aP(:,s-4:s), 1, 2);
        J0 = approxJ(dU, dP);
    end
    
    U1 = U0 - J0 \ Phi(U0);
    
    aU(:, s+1) = U1;
    aP(:, s+1) = Phi(U1);
    ErrJ_SQN(s+1) = norm(J0 - J(U0), 'fro');
    if sqrt(h)*norm(Phi(U1)) < 1e-12 ||...
            sqrt(h)*norm(U0-U1) < 1e-12
        break
    end
    U0 = U1;
end

Nit_SQN = s;
Err_SQN = sqrt(h*sum((aU(:,1:s) - aU(:,s+1)).^2));
Phi_SQN = sqrt(h*sum(aP(:,1:s+1).^2));
ErrJ_SQN = ErrJ_SQN(1:s+1);

disp('Sparse Quasi-Newton:')
disp(s)

%% Sparse Broyden
U0 = ones(N-1, 1);

aU = zeros(N-1, 50+1, 'double');
aP = zeros(N-1, 50+1, 'double');
ErrJ_SB = zeros(1, 50+1, 'double');

aU(:,1) = U0;
aP(:,1) = Phi(U0);
J0 = A(U0);
ErrJ_SB(1) = norm(J0 - J(U0), 'fro');
for s = 1:50
    
    U1 = U0 - J0 \ Phi(U0);
    
    dU = U1 - U0;
    dP = Phi(U1) - Phi(U0);
    J0 = updateJ(J0, dU, dP);
    
    aU(:, s+1) = U1;
    aP(:, s+1) = Phi(U1);
    ErrJ_SB(s+1) = norm(J0 - J(U0), 'fro');
    if sqrt(h)*norm(Phi(U1)) < 1e-12 ||...
            sqrt(h)*norm(U0-U1) < 1e-12
        break
    end
    U0 = U1;
end

Nit_SB = s;
Err_SB = sqrt(h*sum((aU(:,1:s) - aU(:,s+1)).^2));
Phi_SB = sqrt(h*sum(aP(:,1:s+1).^2));
ErrJ_SB = ErrJ_SB(1:s+1);

disp('Sparse Broyden:')
disp(s)

%% plot solution
% figure
% U_exact = arrayfun(ue, x);
% plot(x, U_exact, x, [ue(0); U1; ue(1)])
% legend('exact', 'numerical')

%% compare error
figure
hold on
plot(1:Nit_N, Err_N', '.-')
plot(1:Nit_P, Err_P', '.-')
plot(1:Nit_B, Err_B', '.-')
plot(1:Nit_SP, Err_SP', '.-')
plot(1:Nit_SN, Err_SN', '*-')
plot(1:Nit_SQN, Err_SQN', '*-')
plot(1:Nit_SB, Err_SB', '*-')
set(gca, 'Yscale', 'log')
xlabel('Nit')
ylabel('error(u)')
legend('Newton', 'Picard', 'Broyden', ...
    'Stab-Picard', 'Sparse Newton', 'SQN',...
    'Sparse Broyden')

figure
hold on
plot(0:Nit_N, Phi_N', '.-')
plot(0:Nit_P, Phi_P', '.-')
plot(0:Nit_B, Phi_B', '.-')
plot(0:Nit_SP, Phi_SP', '.-')
plot(0:Nit_SN, Phi_SN', '*-')
plot(0:Nit_SQN, Phi_SQN', '*-')
plot(0:Nit_SB, Phi_SB', '*-')
set(gca, 'Yscale', 'log')
xlabel('Nit')
ylabel('\phi(u)')
legend('Newton', 'Picard', 'Broyden', ...
    'Stab-Picard', 'Sparse Newton', 'SQN',...
    'Sparse Broyden')
