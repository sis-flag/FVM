function J = approxJ(nzInd, dU, dP)

N = size(dU, 1);
deepth = size(dU, 2);

J = spalloc(N, N, 13*N);
for k =1:N
    nzi = nzInd{k};
    if length(nzi) > deepth
        error('I need more')
    end
    
    tA = dU(nzi, :)';
    tb = dP(k, :)';
    J(k, nzi) = (tA \ tb)';
end

end