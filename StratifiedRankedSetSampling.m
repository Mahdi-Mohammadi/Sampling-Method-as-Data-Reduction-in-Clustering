function [ data,varAvg,stdAvg,meanAvg ] = StratifiedRankedSetSampling(dataWithClassName,samplingLoop,stratifyNumber,sampleSize )
dataWithClassName = sort(dataWithClassName);
stratifiedIndex = quantile(dataWithClassName,(1:stratifyNumber-1)/stratifyNumber);
dataMean = zeros(samplingLoop,1);
dataAll = zeros(sampleSize,samplingLoop*2);
selectedData = zeros(sampleSize,2);

for i = 1:samplingLoop
    if (i-1) == 0
            minValue = -Inf;
    else
        minValue = stratifiedIndex(i-1);
    end
    maxValue = stratifiedIndex(i);
    if i==samplingLoop
            maxValue = Inf;
    end
    tempIndex = dataWithClassName(:,1)>=minValue & dataWithClassName(:,1)<maxValue;
    tempDataStratified = dataWithClassName(tempIndex,:);
    for j = 1:sampleSize
        tempDataRanked = sort(tempDataStratified(randsample(1:length(tempDataStratified), sampleSize, true),:));
        selectedData(j,1) = tempDataRanked(j,1);
        selectedData(j,2) = tempDataRanked(j,2);
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