% Write binary STL from face/vertex data
%tmpvol = false(20,20,20); % Empty voxel volume
%tmpvol(8:12,8:12,5:15) = 1; % Turn some voxels on
%fv = isosurface(~tmpvol, 0.5); % Make patch w. faces "out"

tmpvol=padarray(shape_matrix,[1,1,1],0);
fv = isosurface(~tmpvol, 0.5); % Make patch w. faces "out"
stlwrite('test.stl',fv) % Save to binary .stl