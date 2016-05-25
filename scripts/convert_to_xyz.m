xyz_matrix=zeros(sum(shape_matrix(:)),3);
loop_index=1;
for i=1:size(shape_matrix,1)
    for j=1:size(shape_matrix,2)
        for k=1:size(shape_matrix,3)
            
            if shape_matrix(i,j,k)==1;
            xyz_matrix(loop_index,1:3)=[x(i) y(j) z(k)];
            loop_index=loop_index+1;
            end
        end
    end
    i/size(shape_matrix,1)
end
