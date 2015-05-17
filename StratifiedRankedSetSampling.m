function [ data,varAvg,stdAvg,meanAvg ] = StratifiedRankedSetSampling(dataWithClassName,samplingLoop,stratifyNumber,sampleSize,type )
% type of sampling : 
%                   - 1 = Diagonal
%                   - 2 = Middle
%                   - 3 = First point
%                   - 4 = Last point
    
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
        if type == 1 % diagonal
            selectedData(j,1) = tempDataRanked(j,1);
            selectedData(j,2) = tempDataRanked(j,2);
        elseif type == 2 % middle
            middle = round(length(tempDataRanked(:,1))/2);
            selectedData(j,1) = tempDataRanked(middle,1);
            selectedData(j,2) = tempDataRanked(middle,2);
        elseif type == 3 % first point
            selectedData(j,1) = tempDataRanked(1,1);
            selectedData(j,2) = tempDataRanked(1,2);
        elseif type == 4 % last point
            selectedData(j,1) = tempDataRanked(length(tempDataRanked(:,1)),1);
            selectedData(j,2) = tempDataRanked(length(tempDataRanked(:,1)),2);
        end
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