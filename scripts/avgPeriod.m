% This script will characterize the average dominant period of a 2D matrix
% by taking the fourier transform of each row and looking for the largest
% frequency component. dimToS

function [avg,std_dev]=avgPeriod(A,realSpaceLength,dimToSlice)

% transpose matrix if period lies along 1st dimension
if dimToSlice==2
    A=A';
end

%% Take FFT of each line in image,look for peak
startInclude=2;
endInclude=ceil(size(A,2)/2);
results=zeros(size(A,1),2);
smoothing_span=.0022*(endInclude-startInclude+1);
for i=1:size(A,1)
    Afft=fft(A(i,:));
    [results(i,1), results(i,2)]=max(smooth(abs(Afft(startInclude:endInclude)),smoothing_span));
end

    %figure,subplot(1,2,1),plot(A(i,:)),subplot(1,2,2),plot(abs(Afft));
    adjusted_avg=(mean(results(:,2))+startInclude-1-1);
    
    avg=realSpaceLength/adjusted_avg;
    std_dev=abs(avg-realSpaceLength/(adjusted_avg-std(results(:,2))));
end