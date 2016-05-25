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
monolayers_per_iteration=1.5*(10/setup.mesh_size);
monolayers_removed_per_iteration=setup.etch_fraction*monolayers_per_iteration;
light_intensity=setup.light_intensity*1E-3; %mW/cm^2*.001W/mw=W/cm^2
smoothing_exponent=2;
depletion_width=setup.depletion_width;
inactive_substrate=1;
dist_fcn_exp=setup.dist_fcn_exp;
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
    
    cd ..
    nCores=dlmread('nCores.txt');
    cd(script_directory)
    
    if ismac==0 && isunix==1
        run_FDTD_script=strcat(['/opt/lumerical/fdtd/bin/fdtd-solutions -nw -run ',strrep(script_directory,' ','\ '),'/lmg_update.lsf']);
        run_FDTD_simulation=strcat(['sh /opt/lumerical/fdtd/bin/fdtd-run-local.sh -n ',num2str(nCores),' ',strrep(script_directory,' ','\ '),'/']);
        extract_FDTD_script=strcat(['/opt/lumerical/fdtd/bin/fdtd-solutions -nw -run ',strrep(script_directory,' ','\ '),'/lmg_extract.lsf']);
        
    else if ismac==1;
            run_FDTD_script=strcat(['/applications/Lumerical/FDTD\ Solutions/FDTD\ Solutions.app/Contents/MacOS/fdtd-solutions -nw -run ',strrep(script_directory,' ', '\ '),'/lmg_update.lsf']);
            run_FDTD_simulation=strcat(['sh /Applications/Lumerical/FDTD\ Solutions/Examples/fdtd-run-local.sh -n ',num2str(nCores),' ',strrep(script_directory,' ','\ '),'/']);
            extract_FDTD_script=strcat(['/applications/Lumerical/FDTD\ Solutions/FDTD\ Solutions.app/Contents/MacOS/fdtd-solutions -nw -run ',strrep(script_directory,' ', '\ '),'/lmg_extract.lsf']);
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
        if setup.source_coherence==0
            for i=1:setup.source_number
                system([run_FDTD_simulation,'temp',num2str(i),'.fsp'])
            end
        elseif setup.source_coherence==1
            system([run_FDTD_simulation,'Iteration\ 0','.fsp'])
        end
        system(extract_FDTD_script);
        
        input_iteration_number=0;
    end
    
    load(strcat('iteration',num2str(input_iteration_number),'.mat'));
    
    % End Program if Last Iteration
    if input_iteration_number==total_iterations,close all,break;end
    
    %% Manipulate Input Iteration Data
    clear pabs
    x=pabs1.x;
    y=pabs1.y;
    z=pabs1.z;
    
    for i=1:setup.source_number
        pabs(i,1:size(x,1),1:size(y,1),1:size(z,1),1,1:size(pabs1.lambda,1),:,:)...
            =reshape(eval(['pabs' num2str(i) '.Pabs']),size(x,1),size(y,1),size(z,1),1,size(pabs1.lambda,1));
    end
    
    %Check to make sure power absorbed is positive
    if min(pabs(:))/max(pabs(:))<-.001,display('The power absorbed was negative in one or more locations. Somethings wrong with the simulation!');
        break;
    end
    
    
    %% Mass Addition Inputs
    if setup.dimension==3
        power_in=light_intensity*(max(x)-min(x))*(max(z)-min(z))*1E4; %w/cm^2*m^2(100cm/m)^2=w
    else
        power_in=light_intensity*(max(x)-min(x))*1E2; %w/cm^2*m(100cm/m)=w/cm
    end
    %Calculate generation rate
    
    generation_rate=zeros(0);
    for i=1:setup.source_number
        for j=1:size(pabs1.lambda)
            if setup.dimension==3
                generation_rate(i,:,:,:,:,j)=power_in*eval(['setup.intensity' num2str(i)])*pabs(i,:,:,:,:,j)/(hc/eval(['pabs' num2str(i) '.lambda(j)']))*1e-6; %j/s*(m^-3)*j*(1e-2 m/cm)^3=cm^-3
            else
                generation_rate(i,:,:,:,:,j)=power_in*eval(['setup.intensity' num2str(i)])*pabs(i,:,:,:,:,j)/(hc/eval(['pabs' num2str(i) '.lambda(j)']))*1e-4; %j/(s*cm)*(m^-2)*j*(1e-2 m/cm)^2=cm^-3
            end
        end
    end
    
    %Weight generation rate by source profile
    if size(pabs1.lambda,1)>1
        source1_profile_normalized=source1_profile/sum(source1_profile(:,2));
        for j=1:size(pabs1.lambda)
            if pabs1.lambda(j)*1E9-source1_profile(j,1)>.01,['Warning, absorption wavelength does not match source1 profile.'...
                    'Check source profile data.'],end
            generation_rate(1,:,:,:,:,j)=generation_rate(1,:,:,:,:,j)*source1_profile_normalized(j,2);
        end
    end
    
    if setup.source_number>1
        if size(pabs2.lambda,1)>1
            source2_profile_normalized=source2_profile/sum(source2_profile(:,2));
            for j=1:size(pabs2.lambda)
                if pabs2.lambda(j)*1E9-source2_profile(j,1)>.01,['Warning, absorption wavelength does not match source2 profile.'...
                        'Check source profile data.'],end
                generation_rate(2,:,:,:,:,j)=generation_rate(2,:,:,:,:,j)*source2_profile_normalized(j,2);
            end
        end
    end
    
    %quick fix for same wavelength sources. Source weighting is included in
    %setup.lsf. Error possible if using two broad band sources with same center
    %wavelength.
    gen_flag=0;
    if setup.source_number==2 && setup.source_coherence==1
        if setup.wavelength1==setup.wavelength2
            gen=generation_rate(1,:,:,:,:,:)/setup.intensity1;
            gen=squeeze(sum(gen,6));
            gen_flag=1;
        end
    end
    
    if gen_flag==0;
        gen=sum(generation_rate(:,:,:,:,:,:),6);
        gen=squeeze(sum(gen,1));
    end
    clear gen_flag
    
    % Mass/iteration
    if setup.dimension==3
        cells_per_iteration=(size(x,1)-1)*(size(z,1)-1)*monolayers_per_iteration;
    else
        cells_per_iteration=(size(x,1)-1)*(size(z,1))*monolayers_per_iteration;
    end
    
    % Shape Matrix
    if input_iteration_number==0
        shape_matrix=int8(imag(index)>0);%squeeze(pabs(1,:,:,:))~=0;
    else
        load(strcat('struct',num2str(input_iteration_number),'.mat'),'shape_matrix');
        shape_matrix=int8(padarray(shape_matrix,[0,size(y,1)-size(shape_matrix,2),0],'post'));
    end
    cell_matrix=int8(shape_matrix+circshift(shape_matrix,[-1 0 0])+circshift(shape_matrix,[0 -1 0])+circshift(shape_matrix,[0 0 -1])+circshift(shape_matrix,[0 -1 -1])+circshift(shape_matrix,[-1 0 -1])+circshift(shape_matrix,[-1 -1 0])+circshift(shape_matrix,[-1 -1 -1])==8);
    
    if setup.dimension==3
        cell_matrix=cell_matrix(1:size(cell_matrix,1)-1,1:size(cell_matrix,2)-1,1:size(cell_matrix,3)-1);
    else
        cell_matrix=cell_matrix(1:size(cell_matrix,1)-1,1:size(cell_matrix,2)-1,1);
    end
    % save substrate cell matrix for etch loop
    if input_iteration_number==0,substrate_cell_matrix=cell_matrix; save('substrate_cell_matrix','substrate_cell_matrix');clear substrate_cell_matrix;end
    
    % Save Figure of Mass Distribution
    % figure
    % surf(x*1E9,z*1E9,y(squeeze(sum(shape_matrix,2)))'*1E9,'linestyle', 'none'),colorbar
    % savefig(strcat('iteration',num2str(input_iteration_number)));
    
    %%  Identify Coordinates within depletion width
    if depletion_width>0;
        flag_depletion_width=0;
        linear_dist=1;
        coord_within_depletion_width=zeros(0);
        mesh_size=setup.mesh_size;
        
        for ind1=-ceil(depletion_width/mesh_size):ceil(depletion_width/mesh_size)
            for ind2=-ceil(depletion_width/mesh_size):ceil(depletion_width/mesh_size)
                for ind3=-ceil(depletion_width/mesh_size):ceil(depletion_width/mesh_size)
                    if sqrt((ind1*mesh_size)^2+(ind2*mesh_size)^2+(ind3*mesh_size)^2)<=depletion_width
                        coord_within_depletion_width=[coord_within_depletion_width; ind1 ind2 ind3];
                    end
                end
            end
        end
        
        % Consolidate coordinates in z dimension for 2D simulations
        if setup.dimension==2
            [coord_within_depletion_width temp2 temp3]=unique(coord_within_depletion_width(:,1:2),'rows');
            identical_z_within_depletion_width=accumarray(temp3,1);
        end
        
    end
    
    %% Mass Addition Loop
    %gen_edited is the average generation rate of each cell
    gen_edited=single(gen+circshift(gen,[-1 0 0])+circshift(gen,[0 -1 0])+circshift(gen,[0 0 -1])+circshift(gen,[-1 -1 0])+circshift(gen,[-1 0 -1])+circshift(gen,[-1 -1 0])+circshift(gen,[-1 -1 -1]))/8;
    if setup.dimension==3
        gen_edited=gen_edited(1:size(gen_edited,1)-1,1:size(gen_edited,2)-1,1:size(gen_edited,3)-1);
    else
        gen_edited=gen_edited(1:size(gen_edited,1)-1,1:size(gen_edited,2)-1,:);
    end
    
    total_cells_added=0;
    cell_addition_loops=0;
    if input_iteration_number==0
        norm_factor=5;
    else
        norm_factor=1;
    end
    
    while total_cells_added < cells_per_iteration
        
        % Create Neighbor Matrix
        
        neighbor_matrix_1=uint8(circshift(cell_matrix,[-1 0 0])+circshift(cell_matrix,[1 0 0])+circshift(cell_matrix,[0 -1 0])+circshift(cell_matrix,[0 1 0])+circshift(cell_matrix,[0 0 -1])+circshift(cell_matrix,[0 0 1]));
        neighbor_matrix_2=uint8(circshift(cell_matrix,[-1 -1 0])+circshift(cell_matrix,[-1 1 0])+circshift(cell_matrix,[1 -1 0])+circshift(cell_matrix,[1 1 0])...
            +circshift(cell_matrix,[1 0 -1])+circshift(cell_matrix,[1 0 1])+circshift(cell_matrix,[-1 0 -1])+circshift(cell_matrix,[-1 0 1])...
            +circshift(cell_matrix,[0 1 -1])+circshift(cell_matrix,[0 1 1])+circshift(cell_matrix,[0 -1 -1])+circshift(cell_matrix,[0 -1 1]));
        neighbor_matrix_3=uint8(circshift(cell_matrix,[1 1 1])+circshift(cell_matrix,[1 1 -1])+circshift(cell_matrix,[1 -1 1])+circshift(cell_matrix,[1 -1 -1])+circshift(cell_matrix,[-1 1 1])+circshift(cell_matrix,[-1 1 -1])+circshift(cell_matrix,[-1 -1 1])+circshift(cell_matrix,[-1 -1 -1]));
        
        
        %imagesc(x*1e9,-y*1e9,neighbor_matrix_1'),colormap(gray),axis image, axis off
        
        % Define Surface
        
        surface_matrix=cell_matrix+int8(2*(0<neighbor_matrix_1 & neighbor_matrix_1<=6 & cell_matrix==0));
        surface_matrix(:,size(surface_matrix,2),:)=0; %% no PB's in Y direction
        
        
        % imagesc(x*1e9,-y*1e9,surface_matrix(:,:,1)'),colormap(gray),axis image, axis off, colorbar
        
        % Calculate average Power Absorbed
        
        gen_avg=(circshift(gen_edited,[1 0 0])+circshift(gen_edited,[-1 0 0])+circshift(gen_edited,[0 -1 0])+circshift(gen_edited,[0 1 0])+circshift(gen_edited,[0 0 1])+circshift(gen_edited,[0 0 -1]))./(single(neighbor_matrix_1)+single(0==neighbor_matrix_1));
        
        %% Calculate gen collected at surface
        if depletion_width>0
            interface_points=surface_matrix==2;
            interface_points_within_depletion_width=zeros(size(interface_points));
            if setup.dimension==3,identical_z_within_depletion_width=ones(size(coord_within_depletion_width,1),1);end
            
            for ind1=1:size(coord_within_depletion_width,1)
                interface_points_within_depletion_width=interface_points_within_depletion_width+...
                    identical_z_within_depletion_width(ind1)*circshift(interface_points,coord_within_depletion_width(ind1,:));
            end
            interface_points_within_depletion_width=(cell_matrix==1).*(neighbor_matrix_1==6).*interface_points_within_depletion_width;
            
            gen_within_depletion_width=gen_edited;
            gen_not_including_suface=gen_edited.*(surface_matrix==1);
            
            for ind1=1:size(coord_within_depletion_width,1)
                
                gen_within_depletion_width=gen_within_depletion_width+...
                    identical_z_within_depletion_width(ind1)*(circshift(gen_not_including_suface./(interface_points_within_depletion_width+(interface_points_within_depletion_width==0)),coord_within_depletion_width(ind1,:)));
                
            end
            gen_within_depletion_width=gen_within_depletion_width.*(neighbor_matrix_1<6);
            % Calculate Probability of Adding Mass if depletion width>0
            cell_addition_probability=(double(gen_within_depletion_width).^dist_fcn_exp)...
                .*(1-(6-double(neighbor_matrix_1))/6)/1^(1/2*smoothing_exponent)...
                .*(1-(12-double(neighbor_matrix_2))/12)/2^(1/2*smoothing_exponent)...
                .*(1-(8-double(neighbor_matrix_3))/8)/3^(1/2*smoothing_exponent)...
                .*(surface_matrix==2);
        end
        
        % Calculate Probability of Adding Mass if depletion width==0
        if depletion_width==0
            cell_addition_probability=(double(gen_avg).^dist_fcn_exp)...
                .*(1-(6-double(neighbor_matrix_1))/6)/1^(1/2*smoothing_exponent)...
                .*(1-(12-double(neighbor_matrix_2))/12)/2^(1/2*smoothing_exponent)...
                .*(1-(8-double(neighbor_matrix_3))/8)/3^(1/2*smoothing_exponent)...
                .*(surface_matrix==2);
        end
        
        
        %Make probability equal for substrate
        if input_iteration_number==0 && inactive_substrate
            cell_addition_probability=cell_addition_probability>0;
        end
        
        cell_addition_probability=double(cell_addition_probability)/max(cell_addition_probability(:))/norm_factor;
        
        % H=surf(mass_addition_probability(:,:,1)','edgecolor','none');
        % set(H, 'linestyle', 'none');
        % set(gcf,'Renderer','Zbuffer');
        % colorbar
        % view([0,90]);
        
        % Evolve Surface
        cells_added=cell_addition_probability>rand(size(surface_matrix));
        cell_matrix=cell_matrix+int8(cells_added);
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
    if setup.dimension==3
        cells_removed_per_iteration=(size(x,1)-1)*(size(z,1)-1)*monolayers_removed_per_iteration;
    else
        cells_removed_per_iteration=(size(x,1)-1)*(size(z,1))*monolayers_removed_per_iteration;
    end
    
    total_cells_removed=0;
    cell_removal_loops=0;
    norm_factor=10;
    full_cell_matrix=cell_matrix;
    load('substrate_cell_matrix');
    cell_matrix=cell_matrix-substrate_cell_matrix;
    
    while total_cells_removed < cells_removed_per_iteration
        
        % Create Neighbor Matrix
        neighbor_matrix_1=uint8(circshift(cell_matrix,[-1 0 0])+circshift(cell_matrix,[1 0 0])+circshift(cell_matrix,[0 -1 0])+circshift(cell_matrix,[0 1 0])+circshift(cell_matrix,[0 0 -1])+circshift(cell_matrix,[0 0 1]));
        neighbor_matrix_2=uint8(circshift(cell_matrix,[-1 -1 0])+circshift(cell_matrix,[-1 1 0])+circshift(cell_matrix,[1 -1 0])+circshift(cell_matrix,[1 1 0])+circshift(cell_matrix,[1 0 -1])+circshift(cell_matrix,[0 1 -1])+circshift(cell_matrix,[-1 0 -1])+circshift(cell_matrix,[0 -1 -1]));
        neighbor_matrix_3=uint8(circshift(cell_matrix,[1 1 1])+circshift(cell_matrix,[1 1 -1])+circshift(cell_matrix,[1 -1 1])+circshift(cell_matrix,[1 -1 -1])+circshift(cell_matrix,[-1 1 1])+circshift(cell_matrix,[-1 1 -1])+circshift(cell_matrix,[-1 -1 1])+circshift(cell_matrix,[-1 -1 -1]));
        
        % Define Surface
        surface_matrix=uint8(cell_matrix+uint8(neighbor_matrix_1<6 & cell_matrix==1));
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
        cell_removal_loops=cell_removal_loops+1;
        
    end
    sum(sum(sum(cell_matrix)))
    %% Write output iteration Data to .mat File
    cell_matrix=cell_matrix+substrate_cell_matrix;
    temp=padarray(cell_matrix,[1,1,1],'post');
    shape_matrix=temp...
        +circshift(temp,[1,0,0])...
        +circshift(temp,[0,1,0])...
        +circshift(temp,[0,0,1])...
        +circshift(temp,[1,1,0])...
        +circshift(temp,[1,0,1])...
        +circshift(temp,[0,1,1])...
        +circshift(temp,[1,1,1]);
    clear temp
    shape_matrix=double(shape_matrix~=0);
    
    
    lower_boundary=1;
    matrix_size=size(shape_matrix);
    
    
    test=permute(shape_matrix,[1,3,2]);
    [trash1, trash2, upper_boundary]=ind2sub(size(test),find(test(:),1,'last'));
    upper_boundary=upper_boundary+1;
    
    shape_matrix=shape_matrix(:,lower_boundary:upper_boundary,:);
    
    y=y(lower_boundary:upper_boundary);
    
    % write structure
    
    temp_name=['struct' num2str(input_iteration_number+1)];
    if setup.dimension==3
        save(temp_name,'shape_matrix','x','y','z');
    else
        save(temp_name,'shape_matrix','x','y');
    end
    clear test shape_matrix
    %% Run FDTD Script Unless Last Iteration
    if input_iteration_number < total_iterations
        system(run_FDTD_script);
        if setup.source_coherence==0
            for i=1:setup.source_number
                system([run_FDTD_simulation,'temp',num2str(i),'.fsp'])
            end
        elseif setup.source_coherence==1
            system([run_FDTD_simulation,'Iteration\ ',num2str(input_iteration_number+1),'.fsp'])
        end
        system(extract_FDTD_script);
    end
    close all
    
end
