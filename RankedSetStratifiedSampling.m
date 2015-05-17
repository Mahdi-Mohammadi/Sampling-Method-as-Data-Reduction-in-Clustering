function [ data,varAvg,stdAvg,meanAvg ] = RankedSetStratifiedSampling(dataWithClassName,samplingLoop,initPopulation,sampleSize )
% This function generate Sample population from data input with Ranked Set
% Stratified Sampling method
% Input:
%   - dataWithClassName : input data with class name
%   - samplingLoop : number of loop for each sampling
%   - initPopulation : number of input data
%   - sampleSize : number of sample size for each sampling
% 
% Return:
%   - matrix of sample data: 
%                               column 1 is sample data generated
%                               column 2 is class name
% 
% Example:
%   RankedSetStratifiedSampling(randi([-10 10],10000,2), 10, 10000,5);
%   - generate 10*5 sample list from input data(randi) 
dataMean = zeros(samplingLoop,1);
selectedData = zeros(sampleSize,2);
dataAll = zeros(sampleSize,samplingLoop*2);
for i = 1:samplingLoop
    for j = 1:sampleSize
        
        
        
        tempData = sort(dataWithClassName(randsample(1:initPopulation, sampleSize, true),:));
        
        
        selectedData(j,1) = tempData(j,1);
        selectedData(j,2) = tempData(j,2);
    end
    dataMean(i) = mean(selectedData(:,1));
    dataAll(1:sampleSize,i) = selectedData(:,1);
    dataAll(1:sampleSize,samplingLoop+i) = selectedData(:,2);
    
end
varAvg =  var(dataMean);
stdAvg = std(dataMean);
meanAvg = mean(dataMean);
data = reshape(dataAll,[],2);
end