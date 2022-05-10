function J1 = updateJ(J0, dU, dP)

N = length(dU)+1;
dJ = zeros(N-1);

tP = dP - J0 * dU;
for k =1:N-1
    asU = sum(dU(max(k-1,1):min(k+1,N-1)).^2);
    dJ(k,k) = tP(k) * dU(k) / asU;
    if k < N-1
        dJ(k,k+1) = tP(k) * dU(k+1) / asU;
    end
    if k > 1
        dJ(k,k-1) = tP(k) * dU(k-1) / asU;
    end
end
J1 = J0 + dJ;

% disp(norm(J0, 'fro'))
% disp(norm(dJ, 'fro'))
end