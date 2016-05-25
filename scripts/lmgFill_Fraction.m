clear all
load setup_parameters
threshold=.5;
image_directory=dir('images*')
cd(image_directory.name)
Files=dir('*iteration_30_height_map.png');
for i=1:size(Files)
    ['processing: ' num2str(i) ' of ' num2str(size(Files,1))]
temp = regexp(Files(i).name,'(\d*)_iteration','tokens');
sim_number(i)=str2double(char(temp{1}));
temp=imread(Files(i).name);
temp=im2bw(temp,threshold);

fill_fraction(i)=sum(temp(:))/size(temp(:),1);



end

cd ..

%% combine like data
number_identical=6;
data=[parameters.intensity1(sim_number)',fill_fraction'];
data=sortrows(data,1);

for i=1:size(data,1)/6
    index=i*6;
    result.intensity(i)=data(index,1);
    result.data(i,:)=[mean(data(index-5:index,2)) std(data(index-5:index,2))];
 
end
%%
errorbar(result.intensity',result.data(:,1),result.data(:,2))
 

csvwrite(['Fill_Fraction' '_' 'threshold=' num2str(threshold) '.csv'],[result.intensity',result.data])