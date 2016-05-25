 load('setup_parameters.mat')
 filename='lmgresults_filter=0.2.mat'
 load(filename)
 
 % combine data from like datasets if necessary
 if 1
 number_identical=6;
 temp1=parameters.intensity1; parameters.intensity1=[];
 temp2=results.period; results.period=[];
 for i=1:number_identical:size(parameters.dimension,2)
       parameters.intensity1=cat(1,parameters.intensity1,temp1(i));
       results.period=cat(2,results.period, [mean(temp2(1,i:i+number_identical-1)); std(temp2(1,i:i+number_identical-1))]);
       
 end
 
 end
 
%scatter(parameters.intensity1,results.period(1,:)*1e9);
errorbar(parameters.intensity1,results.period(1,:)*1e9,results.period(2,:)*1e9,'rx')
xlabel('Source Composition')
ylabel('Period (nm)')
title(['Simulated Period ' num2str(parameters.wavelength1(1))...
    ' + ' num2str(parameters.wavelength2(1)) ' nm sources'])
set(gca,'FontName','Cambria');
set(gca,'FontSize', 20);
axis([0 1 0 max(results.period(1,:))*1e9+50])
savefig(strrep(strrep(filename,'.mat',''),'.',','))
export_fig(strrep(filename,'.',','))
csvwrite(strrep(filename,'.mat','.csv'),cat(2,parameters.intensity1,results.period'*1e9))