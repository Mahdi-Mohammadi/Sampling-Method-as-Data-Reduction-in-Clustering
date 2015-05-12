function [ data,varAvg,stdAvg,meanAvg ] = SimpleRandomSampling(dataWithClassName,samplingLoop,initPopulation,sampleSize )
dataMean = zeros(samplingLoop,1);
dataAll = zeros(sampleSize,samplingLoop*2);
for i = 1:samplingLoop
    data = dataWithClassName(randsample(1:initPopulation, sampleSize, true),:);
    dataMean(i) = mean(data(:,1));
    dataAll(1:sampleSize,i) = data(:,1);
    dataAll(1:sampleSize,samplingLoop+i) = data(:,2);
    
end
varAvg =  var(dataMean);
stdAvg = std(dataMean);
meanAvg = mean(dataMean);
data = reshape(dataAll,[],2);
end