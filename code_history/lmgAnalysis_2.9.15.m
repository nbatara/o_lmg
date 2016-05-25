%% this script will calculate the period of lamellae

load setup_parameters.mat
mkdir('images')
mkdir('figures')
for folder_index=1:sim_number
    cd([num2str(folder_index) '_' sim_name])
    
    %% find and load last .mat file
    matFiles=dir('iteration*.mat');
    iteration_number=0;
    for matFile_index=1:size(matFiles,1)
        temp = regexp(matFiles(matFile_index).name,'iteration(\d*).mat','tokens');
        temp=str2double(temp{1});
        if temp>iteration_number,iteration_number=temp;end
    end
    
    load(['iteration',num2str(iteration_number),'.mat']);
    load(['setup.mat']);

%% Manipulate Input Iteration Data
x=pabs1.x;
y=pabs1.y;
z=pabs1.z;

for i=1:setup.source_number
    pabs(i,1:size(x,1),1:size(y,1),1:size(z,1))=reshape(eval(['pabs' num2str(i) '.Pabs']),size(x,1),size(y,1),size(z,1));
end

%Check to make sure power absorbed is positive    
if min(pabs(:))/max(pabs(:))<-.001,display('The power absorbed was negative in one or more locations. Somethings wrong with the simulation!');
    break;
end

shape_matrix=squeeze(pabs(1,:,:,:))~=0;
%% 
    [X,Y,Z] = meshgrid(x,z,y);
    figure
    test=isosurface(X*1E9,Y*1E9,Z*1E9,permute(shape_matrix,[3 1 2]),0,Z*1E9);
    p=patch(test);
    colormap('default');
    set(p,'EdgeColor','none','FaceColor','interp');
    %isonormals(X*1E9,Y*1E9,Z*1E9,permute(shape_matrix,[3 1 2]),p)

    axis ([min(x) max(x) min(z) max(z) min(y) max(y)]*1E9)
    axis off
    daspect([1,1,1])
    view(2)

    light 
    lighting gouraud
    
    h1=colorbar;
    set(h1,'Color','white','fontsize',20)
    set(gcf,'Color','none') %no background color

%% Calculate Period
for i=1:size(x,1)
    for j=1:size(z,1)
        
        
surface(i,j)=find(shape_matrix(i,:,j),1,'last');
    end
end

[results.period(1,folder_index) results.period(2,folder_index)]=avgPeriod(surface,x(size(x,1)),2)

cd ..
    image_name=[num2str(folder_index) '_iteration_' num2str(iteration_number)];
    savefig(image_name)
    movefile([image_name '.fig'],'figures/');
    export_fig(image_name)
    movefile([image_name '.png'],'images/');
    close
end

save('lmgresults','results');
