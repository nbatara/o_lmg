clear all
load setup_parameters

cd images
Files=dir('*iteration_30_height_map.png');

for i=1:size(Files)
    ['processing: ' num2str(i) ' of ' num2str(size(Files,1))]
temp = regexp(Files(i).name,'(\d*)_iteration','tokens');
sim_number(i)=str2double(char(temp{1}));
temp=imread(Files(i).name);
temp=temp(ceil(size(temp,1)/2),:);
temp=xcorr2(double(temp));

[peaks(i) period(i)]=findpeaks(double(temp(ceil(size(temp,1)/2),ceil(size(temp,2)/2)+1:ceil(size(temp,2)/2)+220)),'NPeaks',2,'MinPeakProminence',max(temp(:))/1000);
%figure,findpeaks(double(temp(ceil(size(temp,1)/2),ceil(size(temp,2)/2)+1:ceil(size(temp,2)/2)+220)),'NPeaks',2,'MinPeakProminence',max(temp(:))/1000);
period(i)=period(i)*parameters.mesh_size(sim_number(i));
end

cd ..

%% combine like data
number_identical=6;
data=[parameters.intensity1(sim_number)',period'];
data=sortrows(data,1);

for i=1:size(data,1)/6
    index=i*6;
    result.intensity(i)=data(index,1);
    result.data(i,:)=[mean(data(index-5:index,2)) std(data(index-5:index,2))];
 
end
%%
errorbar(result.intensity',result.data(:,1),result.data(:,2))


csvwrite('autocoorelation.csv',[result.intensity',result.data])