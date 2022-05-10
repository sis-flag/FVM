function J = approxJ(dU, dP)

N = size(dU,1) + 1;
J = zeros(N-1);

for k =1:N-1
    if k == 1
        tA = dU(k:k+1,:)';
        tb = dP(k,:)';
        J(k,k:k+1) = (tA \ tb)';
    elseif k == N-1
        tA = dU(k-1:k,:)';
        tb = dP(k,:)';
        J(k,k-1:k) = (tA \ tb)';
    else
        tA = dU(k-1:k+1,:)';
        tb = dP(k,:)';
        J(k,k-1:k+1) = (tA \ tb)';
    end
end

end