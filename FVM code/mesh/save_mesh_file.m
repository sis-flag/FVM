function save_mesh_file(file_name, Mesh)

fid = fopen(file_name,'w');

fprintf(fid, 'vertices\n');
fprintf(fid, '          %d\n', Mesh.nP);
for k = 1:Mesh.nP
fprintf(fid, '    %.10f    %.10f\n', Mesh.xp(k), Mesh.yp(k));
end

fprintf(fid, 'cells\n');
fprintf(fid, '          %d\n', length(Mesh.U2P));
for k = 1:length(Mesh.U2P)
    fprintf(fid, '    %d', length(Mesh.U2P{k}));
    for u = Mesh.U2P{k}
        fprintf(fid, '    %d', u);
    end
    fprintf(fid, '\n');
end

fclose(fid);
end