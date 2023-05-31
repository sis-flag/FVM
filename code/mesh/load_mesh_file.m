function Mesh = load_mesh_file(filename)

fileID = fopen(filename, 'r');

str = fgetl(fileID);
len = str2double(fgetl(fileID));
coord = zeros(2, len);
for k = 1:len
str = fgetl(fileID);
temp = textscan(str, '%f');
coord(:,k) = temp{1};
end

str = fgetl(fileID);
len = str2double(fgetl(fileID));
U2P = cell(len, 1);
for k = 1:len
str = fgetl(fileID);
temp = textscan(str, '%f');
U2P{k} = temp{1}(2:end)';
end

fclose(fileID);

Mesh = arrange_polygonal(coord, U2P);
end
