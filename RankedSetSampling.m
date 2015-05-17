function [ data,varAvg,stdAvg,meanAvg ] = RankedSetSampling(dataWithClassName,samplingLoop,initPopulation,sampleSize,type )
% type of sampling : 
%                   - 1 = Diagonal
%                   - 2 = Middle
%                   - 3 = First point
%                   - 4 = Last point
dataMean = zeros(samplingLoop,1);
selectedData = zeros(sampleSize,2);
dataAll = zeros(sampleSize,samplingLoop*2);
for i = 1:samplingLoop
    for j = 1:sampleSize
        tempData = sort(dataWithClassName(randsample(1:initPopulation, sampleSize, true),:));
        if type == 1 % diagonal
            selectedData(j,1) = tempData(j,1);
            selectedData(j,2) = tempData(j,2);
        elseif type == 2 % middle
            middle = round(length(tempData(:,1))/2);
            selectedData(j,1) = tempData(middle,1);
            selectedData(j,2) = tempData(middle,2);
        elseif type == 3 % first point
            selectedData(j,1) = tempData(1,1);
            selectedData(j,2) = tempData(1,2);
        elseif type == 4 % last point
            selectedData(j,1) = tempData(length(tempData(:,1)),1);
            selectedData(j,2) = tempData(length(tempData(:,1)),2);
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