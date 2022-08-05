function n = norm_cell(Mesh, u, p)

if p == inf
    n = max(abs(u));
elseif p > 0
    n = sum(Mesh.area .* abs(u).^p)^(1/p);
else
    n = NaN;
end

end