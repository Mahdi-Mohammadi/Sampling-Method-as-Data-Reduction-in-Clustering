function [ data,varAvg,stdAvg,meanAvg ] = RankedSetSampling(dataWithClassName,samplingLoop,initPopulation,sampleSize )
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