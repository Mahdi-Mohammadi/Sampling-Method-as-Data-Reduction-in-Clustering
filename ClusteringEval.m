function [ result ] = ClusteringEval( dataClassName,clusterDataIndex,numberOfCluster )
% This function Evaluate Clustering
% Input:
%   - dataClassName : True Class Name
%   - clusterDataIndex : Clustering Class Name
%   - numberOfCluster : Number Of Cluster
% 
% Return:
%   - Percent of Error in Clustering
if (numberOfCluster<=6) % If numberOfCluster > 6 then time of evaluation is very long
    P = perms(1:numberOfCluster); % All possible permutations
    errorTemp = zeros(length(P),1);
    for n = 1:length(P);
        error = 0;
        tmpClass = P(n,:);
        for i = 1:length(dataClassName);
            masterClass = dataClassName(i); % Master Class Name
            kmeansClass = clusterDataIndex(i); % Clustering Class Name 
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