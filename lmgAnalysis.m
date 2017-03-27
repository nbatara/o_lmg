%% this script will calculate the period of lamellae
clear all
load setup_parameters.mat
mkdir('images')
mkdir('figures')
addpath('scripts')
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
full_x=pabs1.x;
full_y=pabs1.y;
full_z=pabs1.z;

clear pabs
for i=1:setup.source_number
    pabs(i,1:size(full_x,1),1:size(full_y,1),1:size(full_z,1))=reshape(eval(['pabs' num2str(i) '.Pabs(:,1,1)']),size(full_x,1),size(full_y,1),size(full_z,1));
end

%Check to make sure power absorbed is positive    
if min(pabs(:))/max(pabs(:))<-.001,display('The power absorbed was negative in one or more locations. Somethings wrong with the simulation!');
    break;
end

full_shape_matrix=squeeze(pabs(1,:,:,:))~=0;

%%  Create Image of Cross Section for 2D simulations
if setup.dimension==2
    
    H=surf(full_x*1E9,full_y*1E9,squeeze(+full_shape_matrix(:,:,1))','EdgeColor','none');
    set(H, 'linestyle', 'none');
    set(gcf,'Renderer','Zbuffer');
    colormap('gray');
    caxis([0 1.8]);
    
    %isonormals(X*1E9,Y*1E9,Z*1E9,permute(shape_matrix,[3 1 2]),p)

    axis ([min(full_x) max(full_x) min(full_y) max(full_y)]*1E9)
    axis off
    daspect([1,1,1])
    view(2)
    
    set(gcf,'Color','none') %no background color
end

%%  Create Isosurface of Top Down Structure for 3D simulations
if setup.dimension==3
    [X,Y,Z] = meshgrid(full_x,full_z,full_y);
    figure
    test=isosurface(X*1E9,Y*1E9,Z*1E9,permute(full_shape_matrix,[3 1 2]),0,Z*1E9);
    p=patch(test);
    colormap('default');
    set(p,'EdgeColor','none','FaceColor','interp');
    %isonormals(X*1E9,Y*1E9,Z*1E9,permute(shape_matrix,[3 1 2]),p)

    axis ([min(full_x) max(full_x) min(full_z) max(full_z) min(full_y) max(full_y)]*1E9)
    axis off
    daspect([1,1,1])
    view(2)

    light 
    lighting gouraud
    
    h1=colorbar;
    set(h1,'Color','white','fontsize',20)
    set(gcf,'Color','none') %no background color
end
%% save bw_image of structure height
load(['struct',num2str(iteration_number),'.mat']);
s=size(shape_matrix);
max_index=zeros(s(1),s(3));
for i=1:s(1)
    for j=1:s(3)
[val,max_index(i,j)]=find(shape_matrix(i,:,j),1,'last');
    end
end
min_val=min(max_index(:));
bw_map=max_index-min_val;
bw_map=bw_map/max(bw_map( :));


%% Calculate Period
for i=1:size(full_x,1)
    for j=1:size(full_z,1)
        
        
surface(i,j)=find(shape_matrix(i,:,j),1,'last');
    end
end

% Filter Stunted Lamellae
filter=.2;
surface=surface.*(surface>=(max(surface(:))*filter));

[results.period(1,folder_index) results.period(2,folder_index)]=avgPeriod(surface,full_x(size(full_x,1))-full_x(1),2)

cd ..
    image_name=[num2str(folder_index) '_iteration_' num2str(iteration_number)];
    savefig(image_name)
    movefile([image_name '.fig'],'figures/');
    if setup.dimension==2 %2D
        imwrite(flip(+full_shape_matrix(:,:,1)',1),[image_name '.png']);
    else %3D
        export_fig(image_name)
    end
    movefile([image_name '.png'],'images/');
    imwrite(flip(bw_map',1),[image_name '_height_map.png']);
    movefile([image_name '_height_map.png'],'images/');
    close
end

save(['lmgresults' '_filter=' num2str(filter) '.mat'],'results');


%% Email Notification
myEmail='';
if ismac==1
  system('scutil --get ComputerName > hostname.txt')
    hostname= textread('hostname.txt','%s','whitespace','');
    hostname=strrep(hostname,' ','');
    delete hostname.txt
end

if ismac==0 && isunix==1 
!hostname > hostname.txt
hostname= textread('hostname.txt','%s');
delete hostname.txt
end

setpref('Internet','E_mail',strcat('MATLAB@',hostname))
setpref('Internet','SMTP_Server','mail.caltech.edu')

sendmail(myEmail,'lmgAnalysis Done!',[sprintf('See Directory:\n') strrep(pwd,' ','\ ')])