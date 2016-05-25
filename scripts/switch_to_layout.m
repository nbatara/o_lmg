% This script will identify all .fsp files and switch them  to layout
% mode
% Nicolas Batara 2/26/15
%% Select Top Directory to Convert .fsp files to layout mode
target_directory=uigetdir;
    if ismac==0 && isunix==1
        run_lumerical='/opt/lumerical/fdtd/bin/fdtd-solutions ';
    else if ismac==1;
            run_lumerical='/applications/Lumerical/FDTD\ Solutions/FDTD\ Solutions.app/Contents/MacOS/fdtd-solutions ';
        end
    end
    
    parameters=[' -nw -run ' strrep(pwd,' ', '\ ') '/switch_to_layout.lsf'];
    
%% find all .fsp file in directory

files=rdir([target_directory '/**/*.fsp']);

%% 
for i=1:size(files,1)
    disp(['converting file ' num2str(i) ' of ' num2str(size(files,1))])
    command=[run_lumerical strrep(files(i).name,' ', '\ ') parameters]
    system(command);
end


%% find all .fsp file in directory

files=rdir([target_directory '/**/*_backup.fsp']);

%% remove any backups

for i=1:size(files,1)
    disp(['deleting backup ' num2str(i) ' of ' num2str(size(files,1))])
    command=['rm ' strrep(files(i).name,' ', '\ ')]
    system(command);
end