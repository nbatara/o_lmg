% This script will take user input to
% setup SeTe simulation directories

%% Prompt user for input parameters

clear parameters setup

if strfind(ls,'setup_parameters.mat')>0
    load setup_parameters.mat
else
    
    % Input Dialog Box to Gather User Data, Command Prompt for
    % Non-Graphical Mode Below (Keep Both Consistant!)
    
    if isempty(javachk('desktop'))
        
        options.Resize='on';
        options.WindowStyle='normal';
        options.Interpreter='tex';
        
        defAns={'delete_me' ...    %1
            '1' ...                %2
            '3' ...                %3
            '0' ...                %4
            '0' ...                %5
            '10' ...               %6
            '1000' ...             %7
            '1000' ...             %8
            '1' ...                %9
            '1' ...                %10
            '1' ...                %11
            '10' ...               %12
            '634' ...              %13
            '0' ...                %14
            '1' ...                %15
            '0' ...                %16
            '0'};                  %17
        
        sampleinfo = inputdlg({ ...
            'Simulation Name (no /''s please):'...                  %1
            'Number of Simulations:' ...                            %2
            'Simulaiton Dimension (2 or 3):'...                     %3
            'Depletion Width (nm):' ...                             %4
            'Etch Fraction (0.0-1.0):' ...                          %5
            'Mesh Size (nm):' ...                                   %6
            'Simulation Height (nm):' ...                           %7
            'Simulation Length (nm):'...                            %8
            'Aspect Ratio for 3D Simulation Depth:'...              %9
            'Number of Illumination Sources (1 or 2):'...           %10
            'Temporal Source Coherence? (0=no, 1=yes):'...          %11
            'Total Illumination Power (mW/cm^2):'...                %12
            'Source 1 Wavelength (nm) (0 for profile import:)'...    %13
            'Source 1 Polarization (deg):'...                       %14
            'Relative Intensity of Source 1:'...                    %15
            'Source 2 Wavelength (nm):'...                          %16
            'Source 2 Polarization (deg):'}, ...                    %17
            'Settings', 1,defAns,options);
        
        sim_name=sampleinfo{1};
        sim_number=str2num(sampleinfo{2});
        
        parameters.dimension=str2num(sampleinfo{3});
        if size(parameters.dimension,2)==1, parameters.dimension=parameters.dimension*ones(1,sim_number); end
        
        parameters.depletion_width=str2num(sampleinfo{4});
        if size(parameters.depletion_width,2)==1, parameters.depletion_width=parameters.depletion_width*ones(1,sim_number); end
        
        parameters.etch_fraction=str2num(sampleinfo{5});
        if size(parameters.etch_fraction,2)==1, parameters.etch_fraction=parameters.etch_fraction*ones(1,sim_number); end
        
        parameters.mesh_size=str2num(sampleinfo{6});
        if size(parameters.mesh_size,2)==1, parameters.mesh_size=parameters.mesh_size*ones(1,sim_number); end
        
        parameters.sim_height=str2num(sampleinfo{7});
        if size(parameters.sim_height,2)==1, parameters.sim_height=parameters.sim_height*ones(1,sim_number); end
        
        parameters.sim_length=str2num(sampleinfo{8});
        if size(parameters.sim_length,2)==1, parameters.sim_length=parameters.sim_length*ones(1,sim_number); end
        
        parameters.aspect_ratio=str2num(sampleinfo{9});
        if size(parameters.aspect_ratio,2)==1, parameters.aspect_ratio=parameters.aspect_ratio*ones(1,sim_number); end
        
        parameters.source_number=str2num(sampleinfo{10});
        if size(parameters.source_number,2)==1, parameters.source_number=parameters.source_number*ones(1,sim_number); end
        
        parameters.source_coherence=str2num(sampleinfo{11});
        if size(parameters.source_coherence,2)==1, parameters.source_coherence=parameters.source_coherence*ones(1,sim_number); end
        
        parameters.light_intensity=str2num(sampleinfo{12});
        if size(parameters.light_intensity,2)==1, parameters.light_intensity=parameters.light_intensity*ones(1,sim_number); end
        
        parameters.wavelength1=str2num(sampleinfo{13});
        if size(parameters.wavelength1,2)==1, parameters.wavelength1=parameters.wavelength1*ones(1,sim_number); end
        
        parameters.pol1=str2num(sampleinfo{14});
        if size(parameters.pol1,2)==1, parameters.pol1=parameters.pol1*ones(1,sim_number); end
        
        if str2num(sampleinfo{10})==2;
            
            parameters.intensity1 = str2num(sampleinfo{15});
            if size(parameters.intensity1,2)==1, parameters.intensity1=parameters.intensity1*ones(1,sim_number); end
            
            parameters.wavelength2=str2num(sampleinfo{16});
            if size(parameters.wavelength2,2)==1, parameters.wavelength2=parameters.wavelength2*ones(1,sim_number); end
            parameters.intensity2=1-parameters.intensity1;
            
            parameters.pol2=str2num(sampleinfo{17});
            if size(parameters.pol2,2)==1, parameters.pol2=parameters.pol2*ones(1,sim_number); end
            
        else
            
            parameters.intensity1=ones(1,sim_number);
        end
        
        clear defAns sampleinfo
        
        % Load intensity profile of source1 if wavelength is set to zero
        if parameters.wavelength1==0
            
            flag=0;
            while flag==0
                % Prompt user for file location and load file
                [FileName,PathName] = uigetfile({'*.txt;*.csv'},...
                    'Select source intensity file (wavelength units in nm, odd number of points):',...
                    'MultiSelect','off');
                source1_profile=importdata([PathName FileName]);
                if mod(size(source1_profile,1),2),flag=1; end
            end
            
            
            parameters.wavelength1=median(source1_profile(:,1));
            if size(parameters.wavelength1,2)==1, parameters.wavelength1=parameters.wavelength1(:,ones(1,sim_number)); end
            
            parameters.source_bandwidth1=max(source1_profile(:,1))-min(source1_profile(:,1));
            if size(parameters.source_bandwidth1,2)==1, parameters.source_bandwidth1=parameters.source_bandwidth1*ones(1,sim_number); end
            
            parameters.source_points1=size(source1_profile,1);
            if size(parameters.source_points1,2)==1, parameters.source_points1=parameters.source_points1*ones(1,sim_number); end
            
        else
            parameters.source_bandwidth1=zeros(1,sim_number);
        end
        
        clear defAns sampleinfo
        
    else
        % Command Prompt Setup for Non Graphical Mode
        
        sim_name=input('Enter a name for this simulation set (no /''s please): ','s');
        
        sim_number=input('Enter number of simulations: ');
        
        temp=input('Enter etch fraction (0.0-1.0): ','s');
        parameters.etch_fraction=str2num(temp);
        if size(parameters.etch_fraction,2)==1, parameters.etch_fraction=parameters.etch_fraction*ones(1,sim_number); end
        
        temp=input('Enter Simulaiton Dimension (2 or 3): ','s');
        parameters.dimension=str2num(temp);
        if size(parameters.dimension,2)==1, parameters.dimension=parameters.dimension*ones(1,sim_number); end
        
        temp=input('Enter mesh size (nm): ','s');
        parameters.mesh_size=str2num(temp);
        if size(parameters.mesh_size,2)==1, parameters.mesh_size=parameters.mesh_size*ones(1,sim_number); end
        
        temp=input('Enter simulation height (nm): ','s');
        parameters.sim_height=str2num(temp);
        if size(parameters.sim_height,2)==1, parameters.sim_height=parameters.sim_height*ones(1,sim_number); end
        
        temp=input('Enter simulation length (nm): ','s');
        parameters.sim_length=str2num(temp);
        if size(parameters.sim_length,2)==1, parameters.sim_length=parameters.sim_length*ones(1,sim_number); end
        
        temp=input('Enter aspect ratio for 3D simulation Depth: ','s');
        parameters.aspect_ratio=str2num(temp);
        if size(parameters.aspect_ratio,2)==1, parameters.aspect_ratio=parameters.aspect_ratio*ones(1,sim_number); end
        
        temp=input('Enter number of sources for simulations (1 or 2): ','s');
        parameters.source_number=str2num(temp);
        if size(parameters.source_number,2)==1, parameters.source_number=parameters.source_number*ones(1,sim_number); end
        
        temp=input('Temporal Source Coherence? (0=no, 1=yes): ','s');
        parameters.source_coherence=str2num(temp);
        if size(parameters.source_coherence,2)==1, parameters.source_coherence=parameters.source_coherence*ones(1,sim_number); end
        
        temp=input('Total Illumination Power (mW/cm^2): ','s');
        parameters.light_intensity=str2num(temp);
        if size(parameters.light_intensity,2)==1, parameters.light_intensity=parameters.light_intensity*ones(1,sim_number); end
        
        temp=input(['Enter wavelength of source_1 for simulation(s) (nm): '],'s');
        parameters.wavelength1=str2num(temp);
        if size(parameters.wavelength1,2)==1, parameters.wavelength1=parameters.wavelength1*ones(1,sim_number); end
        
        temp=input(['Enter polarization of source_1 for simulation(s) (nm): '],'s');
        parameters.pol1=str2num(temp);
        if size(parameters.pol1,2)==1, parameters.pol1=parameters.pol1*ones(1,sim_number); end
        
        if parameters.source_number(1)==2; % Relative intensity of two sources
            temp=input(['Enter intensity of source_1 for simulation(s) (0.0 - 1.0): '],'s');
            parameters.intensity1 = str2num(temp);
            if size(parameters.intensity1,2)==1, parameters.intensity1=parameters.intensity1*ones(1,sim_number); end
            
            temp=input(['Enter wavelength of source2 for simulation(s) (nm): '],'s');
            parameters.wavelength2=str2num(temp);
            if size(parameters.wavelength2,2)==1, parameters.wavelength2=parameters.wavelength2*ones(1,sim_number); end
            
            temp=input(['Enter polarization of source2 for simulation(s) (deg): '],'s');
            parameters.pol2=str2num(temp);
            if size(parameters.pol2,2)==1, parameters.pol2=parameters.pol2*ones(1,sim_number); end
            
            parameters.intensity2=1-parameters.intensity1;
            
        else
            parameters.intensity1=ones(1,sim_number);
        end
        
    end
    
end

parameters
%% Setup Simulation Directories
starting_dir=pwd;
dir_name=[sim_name '_' datestr(datetime,'mm.dd.yy,HH.MM.SS')];
mkdir(dir_name)
cd(dir_name)
save('setup_parameters','parameters','sim_name','sim_number')

%transpose arrays within parameters structure and write to .txt file
field_name=fieldnames(parameters);
for j=1:size(field_name,1)
    eval(['temp.' field_name{j} '=' 'parameters.' field_name{j} ''';']);
end
export(struct2dataset(temp),'file','setup_parameters.txt')
clear temp

% save source1_profile if using
if isfield(parameters,'source_bandwidth1')
    if parameters.source_bandwidth1~= 0,save('source1_profile','source1_profile'),end
end
copyfile([starting_dir '/run_lmg.m'],pwd)
copyfile([starting_dir '/lmgAnalysis.m'],pwd)
%Copy files to subdirectories
for i=1:sim_number
    mkdir([num2str(i) '_' sim_name])
    cd([num2str(i) '_' sim_name])
    %Copy files
    copyfile([starting_dir '/lmg.m'],pwd)
    copyfile([starting_dir '/lmg.lsf'],pwd)
    copyfile([starting_dir '/setup.lsf'],pwd)
    copyfile([starting_dir '/lmg.fsp'],pwd)
    %Save simulation parameters
    
    for j=1:size(field_name,1)
        eval(['setup.' field_name{j} '=' 'parameters.' field_name{j} '(i);']);
    end
    save('setup','setup');
    
    % Save source1_profile if using
    if parameters.source_bandwidth1~=0
        temp=source1_profile;
        
        if size(source1_profile,2)-1==sim_number
            i2=i;
        elseif size(source1_profile,2)-1==1
            i2=1;
        else
            ['Warining: Source Profile size does not match number of simulations']
        end
        source1_profile=[source1_profile(:,1) source1_profile(:,i2+1)];
        
        save('setup','setup','source1_profile');
        source1_profile=temp;
        clear temp
    end
    
    display(['setting up simulation ' num2str(i) ' of ' num2str(sim_number)]);
    system('/Applications/Lumerical/FDTD\ Solutions/FDTD\ Solutions.app/Contents/MacOS/fdtd-solutions -nw -run setup.lsf > /dev/null 2>&1');
    
    
    cd ..
end
cd ..

clear setup;


