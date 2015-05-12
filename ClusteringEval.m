function [ result ] = ClusteringEval( dataClassName,clusterDataIndex,numberOfCluster )
%CLUSTERINGEVAL Summary of this function goes here
%   Detailed explanation goes here
result = 100;
if (numberOfCluster<7)
    P = perms(1:numberOfCluster);
    errorTemp = zeros(length(P),1);
    for n = 1:length(P);
        error = 0;
        tmpClass = P(n,:);
        for i = 1:length(dataClassName);
            masterClass = dataClassName(i);
            kmeansClass = clusterDataIndex(i);
            if(kmeansClass > numberOfCluster )
                error = error +1;
            elseif(masterClass ~= tmpClass(kmeansClass))
                error = error +1;
            end
        end
        errorTemp(n) = (error*100) / length(dataClassName);
    end
    result = min(errorTemp);
end
end