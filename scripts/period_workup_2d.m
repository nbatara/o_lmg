clear all
load('setup_parameters.mat')
cd images
for i=1:size(parameters.dimension,2)
    data=imread([num2str(i) '_iteration_30_height_map.png'])
    dft_data=abs(fft(data(1,:)))
    [val ind]=max(dft_data(2:size(dft_data,2)))
    result(i)=parameters.sim_length(i)/(ind+1)

end
cd ..