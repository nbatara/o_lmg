% 3D Animation for
% Monte Carlo Photo-Electrochemical Deposition Simulation
total_iteration_number=0; %specify number of iterations to visualize. zero means all

%% Find and load latest .mat in local folder
if total_iteration_number==0
i=0;s=1;
while s>0
    s = max(size(which(strcat('iteration',num2str(i),'.mat'))));
    i=i+1;
end
total_iteration_number=i-2
end
input_iteration_number=0;
%% Load Last .mat file to check maximum power absorbed
load(strcat('iteration',num2str(total_iteration_number),'.mat'));
cmax=max(pabs(:))/6*1E-18;
%% Create Animated 3D Plot

if ismac==1
    writerObj = VideoWriter(strcat('iteration',num2str(total_iteration_number),'_isosurface_slow'),'MPEG-4')
%linux doesn't have MPEG4 codec by default
elseif ismac==0 && isunix==1
    writerObj = VideoWriter(strcat('iteration',num2str(total_iteration_number),'_isosurface_animation'),'Motion JPEG AVI')
end
figure('Position', [100, 100, 560*1.8, 420*1.5]);
writerObj.Quality=50;
open(writerObj);

tic
loop_counter=0
for j = 0:240


if input_iteration_number<total_iteration_number
    'first loop'
if loop_counter==floor(240/(total_iteration_number+1))||input_iteration_number==0;
    'second loop'
    load(strcat('iteration',num2str(input_iteration_number),'.mat'));
    input_iteration_number=input_iteration_number+1
    loop_counter=0;
    
    shape_matrix=(imag(index)>0);
    [X,Y,Z] = meshgrid(x,z,y);
    %[isosurf1,isosurf2,isosurf3]=isosurface(X*1E9,Y*1E9,Z*1E9,permute(shape_matrix,[1 3 2]),0);
    %test=isosurface(X*1E9,Y*1E9,Z*1E9,permute(shape_matrix,[1 3 2]),0,Z*1E9);
    %test=isosurface(X*1E9,Y*1E9,Z*1E9,permute(shape_matrix,[3 1 2]),0,permute(pabs*1E-18,[3 1 2]))%Z*1E9);
    test=isosurface(X*1E9,Y*1E9,Z*1E9,permute(shape_matrix,[3 1 2]),0,Z*1E9);

    %test.facevertexcdata=(test.facevertexcdata-min(test.facevertexcdata(:)))/(max(test.facevertexcdata(:))-min(test.facevertexcdata(:)))
    cla
    p=patch(test);
    isonormals(X*1E9,Y*1E9,Z*1E9,permute(shape_matrix,[3 1 2]),p)
    %set(p,'FaceColor',[0.6 0.6 0.6],'EdgeColor','none');
    
    colormap('default');
    set(p,'EdgeColor','none','FaceColor','interp');
    
    colorbar
    
    %caxis([0 cmax])
    axis ([min(x) max(x) min(z) max(z) min(y) max(y)]*1E9)
    set(gca,'fontsize',16,'color',[1 1 1]*.6)
    axis vis3d
    daspect([1,1,1])
    view(3)
    if input_iteration_number==0
        zoom(1.15)
    end
    light 
    lighting gouraud
    xlabel('x')
    ylabel('z')
    zlabel('y')




end
end


loop_counter=loop_counter+1;
    
%savefig(gcf,strcat('iteration',num2str(input_iteration_number),'_isosurface'));
%hgexport(gcf,strcat('iteration',num2str(input_iteration_number),'_isosurface','.jpg'),hgexport('factorystyle'), 'Format', 'jpeg');

     view(-45,60);
     writeVideo(writerObj,getframe(gcf));
% if 0>=j<=120
%     view(-55+1*j/4,60);
%     writeVideo(writerObj,getframe(gcf));
%     j
% end

% if  121>=j<=240
%     view(-25-1*(j-120)/4,60);
%     writeVideo(writerObj,getframe(gcf));
%     j
% end
end

close(writerObj);

toc  