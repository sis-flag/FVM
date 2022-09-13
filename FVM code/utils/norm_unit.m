function n = norm_unit(Mesh, u, p)

if p == inf
    n = max(abs(u));
elseif p > 0
    n = sum(Mesh.area .* abs(u).^p)^(1/p);
else
    error('p is not valid')
end

end