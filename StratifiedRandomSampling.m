function [ data,varAvg,stdAvg,meanAvg ] = StratifiedRandomSampling(dataWithClassName,samplingLoop,stratifyNumber,sampleSize )
dataWithClassName = sort(dataWithClassName);
stratifiedIndex = quantile(dataWithClassName,(1:stratifyNumber-1)/stratifyNumber);
dataMean = zeros(samplingLoop,1);
dataAll = zeros(sampleSize,samplingLoop*2);
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
    tempData = dataWithClassName(tempIndex,:);
    data = tempData(randsample(1:length(tempData), sampleSize, true),:);
    dataMean(i) = mean(data(:,1));
    dataAll(1:sampleSize,i) = data(:,1);
    dataAll(1:sampleSize,samplingLoop+i) = data(:,2);
end
varAvg =  var(dataMean);
stdAvg = std(dataMean);
meanAvg = mean(dataMean);
data = reshape(dataAll,[],2);
end