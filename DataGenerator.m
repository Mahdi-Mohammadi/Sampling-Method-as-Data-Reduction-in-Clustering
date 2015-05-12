function [ data ] = DataGenerator( initPopulation,numberOfCluster,SIGMA,firstMeanStart,meanInterval )
% This function generate 1d data set with mvnrnd method
% Input:
%   - initPopulation : number of all population to initial
%   - numberOfCluster : number of cluster
%   - SIGMA : sigma
%   - firstMeanStart : start point mean for use an mvnrnd
%   - meanInterval : interval between eche mean of cluster
% 
% Return:
%   - matrix of generated data: 
%                               column 1 is random data generated
%                               column 2 is class name
% 
% Example:
%   DataGenerator(40, 4, 1, 0, 1);
%   - generate 40 random number in 4 cluster with sigma 1, start mean 0 and interval 1 between eche mean of cluster 

for i = 1: numberOfCluster
    numPopulationOfThisCluster = initPopulation/numberOfCluster; % calculate number of population of this cluster
    tempMean = firstMeanStart+(meanInterval*(i-1)); % calculate mean of this cluster
    tempData(:,i)=mvnrnd(tempMean,SIGMA,numPopulationOfThisCluster); % generate data
    tempData(:,i+numberOfCluster)=i; % add class name
end
data = reshape(tempData,[],2); %reshape data
end