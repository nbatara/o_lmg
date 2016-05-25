clear all
load setup_parameters
frequency_min=10
frequency_max=150
cd images
Files=dir('*iteration_30_height_map.png');
for i=1:size(Files)
    ['processing: ' num2str(i) ' of ' num2str(size(Files,1))]
temp = regexp(Files(i).name,'(\d*)_iteration','tokens');
sim_number(i)=str2double(char(temp{1}));
temp=imread(Files(i).name);
temp=temp(ceil(size(temp,1)/2),:);
temp=abs(fft(temp));
x=linspace(0,size(temp,2)-1,size(temp,2));
f = fit(x(frequency_min:frequency_max)',temp(frequency_min:frequency_max)','gauss1');

%figure,hold on,plot(x(2:300),f(x(2:300))),plot(x(2:300),temp(2:300))

frequency(i)=f.b1;
temp2=confint(f);
confidenceInt(i)=temp2(2,2)-temp2(1,2);
period(i)=parameters.sim_length(sim_number(i))/frequency(i);



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
 

csvwrite(['FFT_curve_fit_' num2str(frequency_min) '-' num2str(frequency_max) '.csv'],[result.intensity',result.data])