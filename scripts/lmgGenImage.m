% This script will generate images of the lmg iterations.
Files=dir();
for i=1:size(Files)
    if Files(i).isdir & Files(i).name(1)~='.'
        cd(Files(i).name)
    
        
        matFiles=dir('iteration*.mat');
        if size(matFiles,1)>0
           temp=strsplit(Files(i).name,'_');
           image_dir=[temp{1} '_gen'];
           mkdir(image_dir);
           load('setup.mat');
        end
        
        for j=1:size(matFiles)
            load(matFiles(j).name);
            
            %% Calculate Generation Rate
            clear pabs
            x=pabs1.x;y=pabs1.y;z=pabs1.z;
            
            for i=1:setup.source_number
                pabs(i,1:size(x,1),1:size(y,1),1:size(z,1),1,1:size(pabs1.lambda,1),:,:)...
                    =reshape(eval(['pabs' num2str(i) '.Pabs']),size(x,1),size(y,1),size(z,1),1,size(pabs1.lambda,1));
            end
            
            light_intensity=setup.light_intensity*1E-3; %mW/cm^2*.001W/mw=W/cm^2
            hc=6.626E-34*3E8; % j*m
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
                        generation_rate(i,:,:,:,:,j)=power_in*eval(['setup.intensity' num2str(i)])*pabs(i,:,:,:,j)/(hc/eval(['pabs' num2str(i) '.lambda(j)']))*1e-6; %j/s*(m^-3)*j*(1e-2 m/cm)^3=cm^-3
                    else
                        generation_rate(i,:,:,:,:,j)=power_in*eval(['setup.intensity' num2str(i)])*pabs(i,:,:,:,j)/(hc/eval(['pabs' num2str(i) '.lambda(j)']))*1e-4; %j/(s*cm)*(m^-2)*j*(1e-2 m/cm)^2=cm^-3
                    end
                end
            end
            
            %Weight generation rate by source profile
            if size(pabs1.lambda,1)>1
                source1_profile_normalized=source1_profile/sum(source1_profile(:,2));
                for j=1:size(pabs1.lambda)
                    if pabs1.lambda(j)*1E9-source1_profile(j,1)>1E-6,['Warning, absorption wavelength does not match source profile.'...
                            'Check source profile data.'],break;end
                    generation_rate(1,:,:,j)=generation_rate(1,:,:,j)*source1_profile_normalized(j,2);
                end
            end
            
            %quick fix for same wavelength sources. Source weighting is included in
            %setup.lsf.
            gen_flag=0;
            if setup.source_number==2 && setup.source_coherence==1
                if setup.wavelength1==setup.wavelength2
                    gen=squeeze(generation_rate(1,:,:,:,:,:))/setup.intensity1;
                    gen_flag=1;
                end
            end
            
            if gen_flag==0;
                gen=squeeze(sum(generation_rate(:,:,:,:,:,:),1));
            end
            clear gen_flag
            
            
            
            %% Generate and Save Plot
            close all
            FS=26;
            figure('Color','none','position',[1 1 1260 800]);
            surf(x*1e6,y*1e6,squeeze(gen(:,:,1))','EdgeColor','none'),view(2)
            axis tight
            xlabel('microns','Fontsize',FS)
            ylabel('microns','Fontsize',FS)
            colorbar
            set(gca,'TickDir','out','Fontsize',FS)
            
            image_name=[strrep(matFiles(j).name,'.mat','') '_gen'];
            export_fig(image_name)
            movefile([image_name '.png'],image_dir);
        end
        cd ..
    end
end