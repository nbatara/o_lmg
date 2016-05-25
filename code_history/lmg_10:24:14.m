% Monte Carlo Photo-Electrochemical Deposition Simulation
% Coupled With Lumerical's FDTD Software
% Nicolas A. Batara 10/16/14
% Mass Addition Done in Mesh Cell Representaion

%% User Inputs
total_iterations=30;
material=1; % 1:SeTe,2:CdTe 3:a-Se 4:PbSe 5: a-Se:PbSe, 8:1
monolayers_per_iteration=1;
light_intensity=.010; %W/cm^2
smoothing_exponent=1;

load setup.mat
%% Material and Physical Constants

if material==1 % SeTe
delta_t=1; % s
se_molw=78.96; % g/mol
te_molw=127.60; % g/mol
composition= .5; % se(1-x)te(x)
sete_molw= se_molw*(1-composition)+ te_molw*(1-composition);    
n=8; % mol electrons/(mole sete)
density= 4.5*1E6;% g/cm^3*conversion=g/m^3
tao=1E-6; % estimated hole and electron carrier lifetimes
ni=1E10  ; %cm^-3
n0=ni; %cm^-3
p0=ni; %cm^-3
%nk_data=importdata('SeTe_newnkdata.txt'); %import nk data in form [wavelength(nm),n,k]
%material_index=interp1(nk_data(:,1),nk_data(:,2),lambda*10^9,'spline')+interp1(nk_data(:,1),nk_data(:,3),lambda*10^9,'spline')*sqrt(-1);
end
 
if material==2 % CdTe
density= 6.2*1E6;% g/cm^3*conversion=g/m^3
tao=1E-9; % estimated hole and electron carrier lifetimes
ni=1E11; %cm^-3 (Electrodeposited, Miyake Et. Al. J Electroanalytical Chemistry 562 (2004) 247)
n0=ni; %cm^-3
p0=ni; %cm^-3
end

if material==3 % a-Si
density= 4.189E6;% g/m^3
tao=25E-6; % estimated hole and electron carrier lifetimes
ni=1E8; %cm^-3 (Electrodeposited, Miyake Et. Al. J Electroanalytical Chemistry 562 (2004) 247)
n0=ni; %cm^-3
p0=ni; %cm^-3
end

if material==4 %PbSe
density=8.27E6;% g/m^3
tao=1E-6; % estimated hole and electron carrier lifetimes
ni=10^18; %cm^-3 (Epitaxial, Zemel Et. Al. PhysRev 140 (1965) 330)
n0=ni; %cm^-3
p0=ni; %cm^-3
end  

if material==5 %a-Se:PbSe , 8:1
density=8/9*4.189E6+1/9*8.27E6;% g/m^3
tao=8/9*25E-6+1/9*1E-6; % estimated hole and electron carrier lifetimes
ni=8/9*1E8+1/9*10^18; %cm^-3 (Epitaxial, Zemel Et. Al. PhysRev 140 (1965) 330)
n0=ni; %cm^-3
p0=ni; %cm^-3
end  

% Physical Constants
hc=6.626E-34*3E8; % j*m
NA=6.022E23;

%% Number of Iterations Loop
input_iteration_number=-1;
while input_iteration_number <= total_iterations
    
%% Find and load latest .mat in local folder
script_directory=pwd;
if ismac==0 && isunix==1 
run_FDTD_script=strcat(['/opt/lumerical/fdtd/bin/fdtd-solutions -run ',strrep(script_directory,' ','\ '),'/lmg.lsf']);
else if ismac==1;
run_FDTD_script=strcat(['/applications/Lumerical/FDTD\ Solutions/FDTD\ Solutions.app/Contents/MacOS/fdtd-solutions  -run ',strrep(script_directory,' ', '\ '),'/lmg.lsf']);
    end
end
i=0;s=1;
while s>0
    s = max(size(which(strcat('iteration',num2str(i),'.mat'))));
    i=i+1;
end
input_iteration_number=i-2;

% if no .mat file, must run .fsp script to generate
if input_iteration_number==-1
    system(run_FDTD_script);
    input_iteration_number=0;   
end
    
load(strcat('iteration',num2str(input_iteration_number),'.mat'));
if min(pabs(:))/max(pabs(:))<-.001,display('The power absorbed was negative in one or more locations. Somethings wrong with the simulation!');
    break;
end
   
% End Program if Last Iteration
if input_iteration_number==total_iterations,close all,break;end

%% Mass Addition Inputs
power_in=light_intensity*(max(x)-min(x))*(max(z)-min(z))*1E4; %w/cm^2*m^2*(100cm/m)^2
simulation_dimensons=size(pabs);
% Mass/iteration
mass_per_iteration=size(x,1)*size(z,1)*monolayers_per_iteration*(x(2)-x(1))*(y(2)-y(1))*(z(2)-z(1))*density;

% Shape Matrix
shape_matrix=pabs~=0;
cell_matrix=shape_matrix+circshift(shape_matrix,[-1 0 0])+circshift(shape_matrix,[0 -1 0])+circshift(shape_matrix,[0 0 -1])+circshift(shape_matrix,[-1 -1 0])+circshift(shape_matrix,[-1 0 -1])+circshift(shape_matrix,[-1 -1 0])+circshift(shape_matrix,[-1 -1 -1])==8;
cell_matrix=cell_matrix(1:size(cell_matrix,1)-1,1:size(cell_matrix,2)-1,1:size(cell_matrix,3)-1);
figure
surf(x*1E9,z*1E9,y(squeeze(sum(shape_matrix,2)))'*1E9,'linestyle', 'none'),colorbar
savefig(strcat('iteration',num2str(input_iteration_number)));

%% Mass Addition Loop
pabs_edited=pabs+circshift(pabs,[-1 0 0])+circshift(pabs,[0 -1 0])+circshift(pabs,[0 0 -1])+circshift(pabs,[-1 -1 0])+circshift(pabs,[-1 0 -1])+circshift(pabs,[-1 -1 0])+circshift(pabs,[-1 -1 -1])/8;
pabs_edited=pabs_edited(1:size(pabs_edited,1)-1,1:size(pabs_edited,2)-1,1:size(pabs_edited,3)-1);

total_mass_added=0;

while total_mass_added < mass_per_iteration

% Create Neighbor Matrix
tic
neighbor_matrix_1=uint8(circshift(cell_matrix,[-1 0 0])+circshift(cell_matrix,[1 0 0])+circshift(cell_matrix,[0 -1 0])+circshift(cell_matrix,[0 1 0])+circshift(cell_matrix,[0 0 -1])+circshift(cell_matrix,[0 0 1]));
neighbor_matrix_2=uint8(circshift(cell_matrix,[-1 -1 0])+circshift(cell_matrix,[-1 1 0])+circshift(cell_matrix,[1 -1 0])+circshift(cell_matrix,[1 1 0])+circshift(cell_matrix,[1 0 -1])+circshift(cell_matrix,[0 1 -1])+circshift(cell_matrix,[-1 0 -1])+circshift(cell_matrix,[0 -1 -1]));
neighbor_matrix_3=uint8(circshift(cell_matrix,[1 1 1])+circshift(cell_matrix,[1 1 -1])+circshift(cell_matrix,[1 -1 1])+circshift(cell_matrix,[1 -1 -1])+circshift(cell_matrix,[-1 1 1])+circshift(cell_matrix,[-1 1 -1])+circshift(cell_matrix,[-1 -1 1])+circshift(cell_matrix,[-1 -1 -1]));
toc

%imagesc(x*1e9,-y*1e9,neighbor_matrix_1'),colormap(gray),axis image, axis off

% Define Surface
tic
surface_matrix=uint8(cell_matrix+2*(0<neighbor_matrix_1 & neighbor_matrix_1<=6 & cell_matrix==0));
surface_matrix(:,size(surface_matrix,2),:)=0; %% no PB's in Y direction 
toc

% imagesc(x*1e9,-y*1e9,surface_matrix(:,:,1)'),colormap(gray),axis image, axis off, colorbar

% Calculate average Power Absorbed
tic
pabs_avg=(circshift(pabs_edited,[1 0 0])+circshift(pabs_edited,[-1 0 0])+circshift(pabs_edited,[0 -1 0])+circshift(pabs_edited,[0 1 0])+circshift(pabs_edited,[0 0 1])+circshift(pabs_edited,[0 0 -1]))./(double(neighbor_matrix_1)+double(0==neighbor_matrix_1));

toc

% Calculate Probability of Adding Mass
tic
mass_addition_probability=(1+((n0+p0)*power_in*pabs_avg/10^4/(hc/lambda(1))...
    *tao+((power_in*pabs_avg/10^4/(hc/lambda(1))*tao).^2)*ni^-2))...
    .*(1-(6-double(neighbor_matrix_1))/6/1^(1/2*smoothing_exponent))...
    .*(1-(12-double(neighbor_matrix_2))/12/2^(1/2*smoothing_exponent))...
    .*(1-(8-double(neighbor_matrix_3))/8/3^(1/2*smoothing_exponent))...
    .*(surface_matrix==2);
mass_addition_probability=mass_addition_probability/max(mass_addition_probability(:))/10;
toc

% H=surf(mass_addition_probability(:,:,1)','edgecolor','none');
% set(H, 'linestyle', 'none');
% set(gcf,'Renderer','Zbuffer');
% colorbar
% view([0,90]);

%% Evolve Surface

tic
mass_added=mass_addition_probability>rand(size(surface_matrix));
cell_matrix=cell_matrix+mass_added;
pabs_edited=pabs_edited+mass_added.*pabs_avg;
total_mass_added=total_mass_added+sum(mass_added(:))*density*(x(2)-x(1))*(y(2)-y(1))*(z(2)-z(1));
toc

% H=surf(x*1E9,y*1E9,double(mass_added(:,:,1))','edgecolor','none');
% set(H, 'linestyle', 'none');
% set(gcf,'Renderer','Zbuffer');
% colorbar
% view([0,90]);
end

%% Write output iteration Data to Text File
s1=size(shape_matrix,1);s2=size(shape_matrix,2);s3=size(shape_matrix,1);
shape_matrix(1:s1-1,1:s2-1,1:s3-1)=cell_matrix;
shape_matrix(2:s1,1:s2-1,1:s3-1)=cell_matrix;
shape_matrix(1:s1-1,2:s2,1:s3-1)=cell_matrix;
shape_matrix(1:s1-1,1:s2-1,2:s3)=cell_matrix;
shape_matrix(2:s1,2:s2,1:s3-1)=cell_matrix;
shape_matrix(2:s1,1:s2-1,2:s3)=cell_matrix;
shape_matrix(1:s1-1,2:s2,2:s3)=cell_matrix;
shape_matrix(2:s1,2:s2,2:s3)=cell_matrix;

lower_boundary=1;
matrix_size=size(shape_matrix);

tic
test=permute(shape_matrix,[1,3,2]);
[trash1, trash2, upper_boundary]=ind2sub(size(test),find(test(:),1,'last'));
upper_boundary=upper_boundary+1; 
clear test
        

% write header
fid=strcat('iteration',num2str(input_iteration_number+1),'.txt');
fprintf(fopen(fid,'w'),'%s',strcat([num2str(matrix_size(1)),' ',num2str(x(1)),' ',num2str(x(matrix_size(1)))]));
fclose('all');
fprintf(fopen(fid,'a+'),'\n%s',strcat([num2str(upper_boundary-lower_boundary+1),' ',num2str(y(lower_boundary)),' ',num2str(y(upper_boundary))]));
fclose('all');
fprintf(fopen(fid,'a+'),'\n%s\n',strcat([num2str(matrix_size(3)),' ',num2str(z(1)),' ',num2str(z(matrix_size(3)))]));
fclose('all');

% write structure
index_edited_for_write=shape_matrix(:,lower_boundary:upper_boundary,:);
dlmwrite(fid,transpose(real(index_edited_for_write(1:size(index_edited_for_write(:))))),'-append','delimiter',' ');

toc
%% Run FDTD Script Unless Last Iteration
if input_iteration_number < total_iterations
    system(run_FDTD_script);
end
close all

end
