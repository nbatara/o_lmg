% Monte Carlo Photo-Electrochemical Deposition Simulation
% Coupled With Lumerical's FDTD Software
% Nicolas A. Batara 12/11/14
% Mass Addition Done in Mesh Cell Representaion

load setup.mat
if exist('working.mat'),disp('This simulation might already be already running. Ending Early.'),...
        return, else working=1; save('working','working'); end
%% User Inputs
total_iterations=30;
material=1; % 1:SeTe,2:CdTe 3:a-Se 4:PbSe 5: a-Se:PbSe, 8:1
monolayers_per_iteration=1;
monolayers_removed_per_iteration=.75;
light_intensity=.010; %W/cm^2
smoothing_exponent=1;

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
   
% End Program if Last Iteration
if input_iteration_number==total_iterations,close all,break;end

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


%% Mass Addition Inputs
power_in=light_intensity*(max(x)-min(x))*(max(z)-min(z))*1E4; %w/cm^2*m^2(100cm/m)^2=w

%Calculate generation rate

generation_rate=zeros(0);
for i=1:setup.source_number
generation_rate(i,:,:,:)=power_in*eval(['setup.intensity' num2str(i)])*pabs(i,:,:,:)/(hc/eval(['pabs' num2str(i) '.lambda']))*1e-6; %j/s*(m^-3)*j*(1e-2 m/cm)^3=cm^-3
end
gen=squeeze(sum(generation_rate,1));

% Mass/iteration
cells_per_iteration=(size(x,1)-1)*(size(z,1)-1)*monolayers_per_iteration;

% Shape Matrix
shape_matrix=squeeze(pabs(1,:,:,:))~=0;
cell_matrix=shape_matrix+circshift(shape_matrix,[-1 0 0])+circshift(shape_matrix,[0 -1 0])+circshift(shape_matrix,[0 0 -1])+circshift(shape_matrix,[-1 -1 0])+circshift(shape_matrix,[-1 0 -1])+circshift(shape_matrix,[-1 -1 0])+circshift(shape_matrix,[-1 -1 -1])==8;
cell_matrix=cell_matrix(1:size(cell_matrix,1)-1,1:size(cell_matrix,2)-1,1:size(cell_matrix,3)-1);

% Save Figure of Mass Distribution
% figure
% surf(x*1E9,z*1E9,y(squeeze(sum(shape_matrix,2)))'*1E9,'linestyle', 'none'),colorbar
% savefig(strcat('iteration',num2str(input_iteration_number)));

%% Mass Addition Loop
%gen_edited is the average generation rate of each cell
gen_edited=(gen+circshift(gen,[-1 0 0])+circshift(gen,[0 -1 0])+circshift(gen,[0 0 -1])+circshift(gen,[-1 -1 0])+circshift(gen,[-1 0 -1])+circshift(gen,[-1 -1 0])+circshift(gen,[-1 -1 -1]))/8;
gen_edited=gen_edited(1:size(gen_edited,1)-1,1:size(gen_edited,2)-1,1:size(gen_edited,3)-1);

total_cells_added=0;
cell_addition_loops=0;
if input_iteration_number==0
    norm_factor=10;
else
    norm_factor=1;
end

while total_cells_added < cells_per_iteration

% Create Neighbor Matrix

neighbor_matrix_1=uint8(circshift(cell_matrix,[-1 0 0])+circshift(cell_matrix,[1 0 0])+circshift(cell_matrix,[0 -1 0])+circshift(cell_matrix,[0 1 0])+circshift(cell_matrix,[0 0 -1])+circshift(cell_matrix,[0 0 1]));
neighbor_matrix_2=uint8(circshift(cell_matrix,[-1 -1 0])+circshift(cell_matrix,[-1 1 0])+circshift(cell_matrix,[1 -1 0])+circshift(cell_matrix,[1 1 0])+circshift(cell_matrix,[1 0 -1])+circshift(cell_matrix,[0 1 -1])+circshift(cell_matrix,[-1 0 -1])+circshift(cell_matrix,[0 -1 -1]));
neighbor_matrix_3=uint8(circshift(cell_matrix,[1 1 1])+circshift(cell_matrix,[1 1 -1])+circshift(cell_matrix,[1 -1 1])+circshift(cell_matrix,[1 -1 -1])+circshift(cell_matrix,[-1 1 1])+circshift(cell_matrix,[-1 1 -1])+circshift(cell_matrix,[-1 -1 1])+circshift(cell_matrix,[-1 -1 -1]));


%imagesc(x*1e9,-y*1e9,neighbor_matrix_1'),colormap(gray),axis image, axis off

% Define Surface

surface_matrix=uint8(cell_matrix+2*(0<neighbor_matrix_1 & neighbor_matrix_1<=6 & cell_matrix==0));
surface_matrix(:,size(surface_matrix,2),:)=0; %% no PB's in Y direction 


% imagesc(x*1e9,-y*1e9,surface_matrix(:,:,1)'),colormap(gray),axis image, axis off, colorbar

% Calculate average Power Absorbed

gen_avg=(circshift(gen_edited,[1 0 0])+circshift(gen_edited,[-1 0 0])+circshift(gen_edited,[0 -1 0])+circshift(gen_edited,[0 1 0])+circshift(gen_edited,[0 0 1])+circshift(gen_edited,[0 0 -1]))./(double(neighbor_matrix_1)+double(0==neighbor_matrix_1));



% Calculate Probability of Adding Mass

cell_addition_probability=(1+((n0+p0)*power_in*gen_avg...
    *tao+((power_in*gen_avg*tao).^2)*ni^-2))...
    .*(1-(6-double(neighbor_matrix_1))/6/1^(1/2*smoothing_exponent))...
    .*(1-(12-double(neighbor_matrix_2))/12/2^(1/2*smoothing_exponent))...
    .*(1-(8-double(neighbor_matrix_3))/8/3^(1/2*smoothing_exponent))...
    .*(surface_matrix==2);

cell_addition_probability=cell_addition_probability/max(cell_addition_probability(:))/norm_factor;


% H=surf(mass_addition_probability(:,:,1)','edgecolor','none');
% set(H, 'linestyle', 'none');
% set(gcf,'Renderer','Zbuffer');
% colorbar
% view([0,90]);

% Evolve Surface
cells_added=cell_addition_probability>rand(size(surface_matrix));
cell_matrix=cell_matrix+cells_added;
gen_edited=gen_edited+cells_added.*gen_avg;
total_cells_added=total_cells_added+sum(cells_added(:));
disp([num2str(floor(total_cells_added/cells_per_iteration*100)) '% of cells added'])
cell_addition_loops=cell_addition_loops+1;

% H=surf(x*1E9,y*1E9,double(mass_added(:,:,1))','edgecolor','none');
% set(H, 'linestyle', 'none');
% set(gcf,'Renderer','Zbuffer');
% colorbar
% view([0,90]);
end


%% Mass Removal Loop
cells_removed_per_iteration=(size(x,1)-1)*(size(z,1)-1)*monolayers_removed_per_iteration;
total_cells_removed=0;
cell_removal_loops=0;
norm_factor=10;


while total_cells_removed < cells_removed_per_iteration

% Create Neighbor Matrix
neighbor_matrix_1=uint8(circshift(cell_matrix,[-1 0 0])+circshift(cell_matrix,[1 0 0])+circshift(cell_matrix,[0 -1 0])+circshift(cell_matrix,[0 1 0])+circshift(cell_matrix,[0 0 -1])+circshift(cell_matrix,[0 0 1]));
neighbor_matrix_2=uint8(circshift(cell_matrix,[-1 -1 0])+circshift(cell_matrix,[-1 1 0])+circshift(cell_matrix,[1 -1 0])+circshift(cell_matrix,[1 1 0])+circshift(cell_matrix,[1 0 -1])+circshift(cell_matrix,[0 1 -1])+circshift(cell_matrix,[-1 0 -1])+circshift(cell_matrix,[0 -1 -1]));
neighbor_matrix_3=uint8(circshift(cell_matrix,[1 1 1])+circshift(cell_matrix,[1 1 -1])+circshift(cell_matrix,[1 -1 1])+circshift(cell_matrix,[1 -1 -1])+circshift(cell_matrix,[-1 1 1])+circshift(cell_matrix,[-1 1 -1])+circshift(cell_matrix,[-1 -1 1])+circshift(cell_matrix,[-1 -1 -1]));

% Define Surface
surface_matrix=uint8(cell_matrix+(neighbor_matrix_1<6 & cell_matrix==1));
surface_matrix(:,size(surface_matrix,2),:)=0; % no PB's in Y direction 
surface_matrix(:,1,:)=0; % no PB's in Y direction

% imagesc(x*1e9,-y*1e9,surface_matrix(:,:,1)'),colormap(gray),axis image, axis off, colorbar

% Calculate Probability of Removing Mass

cell_removal_probability=...
    (1-(1-(6-double(neighbor_matrix_1))/6/1^(1/2*smoothing_exponent))...
    .*(1-(12-double(neighbor_matrix_2))/12/2^(1/2*smoothing_exponent))...
    .*(1-(8-double(neighbor_matrix_3))/8/3^(1/2*smoothing_exponent)))...
    .*(surface_matrix==2);

cell_removal_probability=cell_removal_probability/max(cell_removal_probability(:))/norm_factor;


% Evolve Surface
cells_removed=cell_removal_probability>rand(size(surface_matrix));
cell_matrix=cell_matrix-cells_removed;
total_cells_removed=total_cells_removed+sum(cells_removed(:));
disp([num2str(floor(total_cells_removed/cells_removed_per_iteration*100)) '% of cells removed'])
cell_addition_loops=cell_addition_loops+1;

end
%% Write output iteration Data to .mat File

s1=size(shape_matrix,1);s2=size(shape_matrix,2);s3=size(shape_matrix,1);
shape_matrix(1:s1-1,1:s2-1,1:s3-1)=cell_matrix;
shape_matrix(2:s1,1:s2-1,1:s3-1)=cell_matrix;
shape_matrix(1:s1-1,2:s2,1:s3-1)=cell_matrix;
shape_matrix(1:s1-1,1:s2-1,2:s3)=cell_matrix;
shape_matrix(2:s1,2:s2,1:s3-1)=cell_matrix;
shape_matrix(2:s1,1:s2-1,2:s3)=cell_matrix;
shape_matrix(1:s1-1,2:s2,2:s3)=cell_matrix;
shape_matrix(2:s1,2:s2,2:s3)=cell_matrix;

shape_matrix=double(shape_matrix);


lower_boundary=1;
matrix_size=size(shape_matrix);


test=permute(shape_matrix,[1,3,2]);
[trash1, trash2, upper_boundary]=ind2sub(size(test),find(test(:),1,'last'));
upper_boundary=upper_boundary+1; 

shape_matrix=shape_matrix(:,lower_boundary:upper_boundary,:);
y=y(lower_boundary:upper_boundary);

% write structure
index_edited_for_write=shape_matrix(:,lower_boundary:upper_boundary,:);
temp_name=['struct' num2str(input_iteration_number+1)];
save(temp_name,'shape_matrix','x','y','z');

%% Run FDTD Script Unless Last Iteration
if input_iteration_number < total_iterations
    system(run_FDTD_script);
end
close all

end
