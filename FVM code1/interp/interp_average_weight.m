function weight = interp_average_weight(Mesh, ~, ~)
% get interplation weights on mesh points by average interplation

weight = cell(Mesh.nP, 1);
for P = 1:Mesh.nP
    if Mesh.isbdp(P)
        weight{P} = [];
    else
        Nu = length(Mesh.P2U{P});
        weight{P} = ones(Nu, 1) / Nu;
    end
end

end