clc
clear all;
%--------------------------------------------------------------------------
% Setting
%-------------------------------------------------------------------------
initPopulation = 100000;% Number of Population
sampleSize = 5;% Size of Sample
samplingLoop= 20;% Number of Loop for Sampling
stratifyNumber = samplingLoop;
firstMeanStart = -4;% Mean
meanInterval = 3;% Interval between of 2 Mean
SIGMA =1;% Sigma
numberOfCluster = 4; % Number Of cluster for init(generate data) and k-means
mainLoop = 10; % Number of Loop for All Program
secendLoop = 50; % Number of Loop for each main loop
% Clustering Algorithm
clusteringFunction = @(X,K) kmeans(X, K, ...
    'EmptyAction','singleton', 'Replicates',5, ...
    'Distance','sqeuclidean','MaxIter',500,...
    'Display','off','Options',statset('UseParallel',0));
showMessage = 1; % On/Off Command Window Message
runSimpleRandomSampling = 1; % Run Simple Random
runRankedSetSampling = 1; % Run Ranked Set
runStratifiedRandomSampling = 1; % Run Stratified Random
runStratifiedRankedSetSampling = 1; % Run Stratified Ranked Set
runAllData = 1; % Run Clustering on All Data
%--------------------------------------------------------------------------
% Start Main
%--------------------------------------------------------------------------
fprintf('Start\n');

%--------------------------------------------------------------------------
% Create Zero Vector For SpeedUp
%--------------------------------------------------------------------------
simpleRandomSampleVar = zeros(mainLoop,secendLoop);
simpleRandomSampleStd = zeros(mainLoop,secendLoop);
simpleRandomSampleMean = zeros(mainLoop,secendLoop);
simpleRandomSampleClusteringSumDistance = zeros(mainLoop,secendLoop);
simpleRandomSampleClusterEval = zeros(mainLoop,secendLoop);

rankedSetSampleVar = zeros(mainLoop,secendLoop);
rankedSetSampleStd = zeros(mainLoop,secendLoop);
rankedSetSampleMean = zeros(mainLoop,secendLoop);
rankedSetSampleClusteringSumDistance = zeros(mainLoop,secendLoop);
rankedSetSampleClusterEval = zeros(mainLoop,secendLoop);


stratifiedRandomSampleVar = zeros(mainLoop,secendLoop);
stratifiedRandomSampleStd = zeros(mainLoop,secendLoop);
stratifiedRandomSampleMean = zeros(mainLoop,secendLoop);
stratifiedRandomSampleClusteringSumDistance = zeros(mainLoop,secendLoop);
stratifiedRandomSampleClusterEval = zeros(mainLoop,secendLoop);

stratifiedRankedSetSampleVar = zeros(mainLoop,secendLoop);
stratifiedRankedSetSampleStd = zeros(mainLoop,secendLoop);
stratifiedRankedSetSampleMean = zeros(mainLoop,secendLoop);
stratifiedRankedSetSampleClusteringSumDistance = zeros(mainLoop,secendLoop);
stratifiedRankedSetSampleClusterEval = zeros(mainLoop,secendLoop);


allDataClusteringSumDistance = zeros(mainLoop,1);
allDataClusterEval = zeros(mainLoop,1);
%--------------------------------------------------------------------------
% Start of Main Loop
%--------------------------------------------------------------------------
for i = 1: mainLoop
    if showMessage
        fprintf(sprintf('Main Loop %d Start \n',i));
    end
    %--------------------------------------------------------------------------
    % Generate Base Data
    %--------------------------------------------------------------------------
    dataWithClassName = DataGenerator( initPopulation,numberOfCluster,SIGMA,firstMeanStart,meanInterval );
    %dataRow = dataWithClassName(:,1);
    
    %--------------------------------------------------------------------------
    % Init Result Matrix Label
    %--------------------------------------------------------------------------
    %resultLabel = {'mainLoop' 'Simple Var' 'Simple Std' 'Simple Mean'};
    
    %--------------------------------------------------------------------------
    % Clustering All Data
    %--------------------------------------------------------------------------
    if runAllData
        %Clustering
        [allDataClusteringIndex,allDataClusteringCenter,allDataClusteringSumDistanceTemp] = clusteringFunction(dataWithClassName(:,1),numberOfCluster);
        allDataClusteringSumDistance(i) = sum(allDataClusteringSumDistanceTemp);
        
        %Clustering Evaluation
        allDataClusterEval(i) = ClusteringEval(dataWithClassName(:,2),allDataClusteringIndex,numberOfCluster);
    end
    
    
    %--------------------------------------------------------------------------
    % Start of Main Loop
    %--------------------------------------------------------------------------
    for j = 1: secendLoop
        if showMessage
            fprintf(sprintf('Secend Loop %d Start \n',j));
        end
        
        %----------------------------------------------------------------------
        %Simple Random Sampling
        %----------------------------------------------------------------------
        if runSimpleRandomSampling
            if showMessage
                %fprintf('Simple Random Sampling\n');
            end
            %Sampling
            [simpleRandomSampleData,simpleRandomSampleVar(i,j),simpleRandomSampleStd(i,j),simpleRandomSampleMean(i,j)] = SimpleRandomSampling(dataWithClassName, samplingLoop, initPopulation, sampleSize);
            
            %Clustering
            [simpleRandomSampleClusteringIndex,simpleRandomSampleClusteringCenter,simpleRandomSampleClusteringSumDistanceTemp] = clusteringFunction(simpleRandomSampleData(:,1),numberOfCluster);
            simpleRandomSampleClusteringSumDistance(i,j) = sum(simpleRandomSampleClusteringSumDistanceTemp);
            
            %Clustering Evaluation
            simpleRandomSampleClusterEval(i,j) = ClusteringEval(simpleRandomSampleData(:,2),simpleRandomSampleClusteringIndex,numberOfCluster);
        end
        
        %----------------------------------------------------------------------
        %Ranked Set Sampling
        %----------------------------------------------------------------------
        if runRankedSetSampling
            if showMessage
                %fprintf('Ranked Set Sampling\n');
            end
            %Sampling
            [rankedSetSampleData,rankedSetSampleVar(i,j),rankedSetSampleStd(i,j),rankedSetSampleMean(i,j)] = RankedSetSampling(dataWithClassName, samplingLoop, initPopulation, sampleSize);
            
            %Clustering
            [rankedSetSampleClusteringIndex,rankedSetSampleClusteringCenter,rankedSetSampleClusteringSumDistanceTemp] = clusteringFunction(rankedSetSampleData(:,1),numberOfCluster);
            rankedSetSampleClusteringSumDistance(i,j) = sum(rankedSetSampleClusteringSumDistanceTemp);
            
            %Clustering Evaluation
            rankedSetSampleClusterEval(i,j) = ClusteringEval(rankedSetSampleData(:,2),rankedSetSampleClusteringIndex,numberOfCluster);
            
        end
        
        %----------------------------------------------------------------------
        %Stratified Random Sampling
        %----------------------------------------------------------------------
        if runStratifiedRandomSampling
            if showMessage
                %fprintf('Stratified Random Sampling\n');
            end
            %Sampling
            [stratifiedRandomSampleData,stratifiedRandomSampleVar(i,j),stratifiedRandomSampleStd(i,j),stratifiedRandomSampleMean(i,j)] = StratifiedRandomSampling(dataWithClassName, samplingLoop, stratifyNumber, sampleSize);
            
            %Clustering
            [stratifiedRandomSampleClusteringIndex,stratifiedRandomSampleClusteringCenter,stratifiedRandomSampleClusteringSumDistanceTemp] = clusteringFunction(stratifiedRandomSampleData(:,1),numberOfCluster);
            stratifiedRandomSampleClusteringSumDistance(i,j) = sum(stratifiedRandomSampleClusteringSumDistanceTemp);
            
            %Run Clustering Evaluation
            stratifiedRandomSampleClusterEval(i,j) = ClusteringEval(stratifiedRandomSampleData(:,2),stratifiedRandomSampleClusteringIndex,numberOfCluster);
            
        end
        
        %----------------------------------------------------------------------
        %Stratified Ranked Set Sampling
        %----------------------------------------------------------------------
        if runStratifiedRankedSetSampling
            if showMessage
                %fprintf('Stratified RankedSet Sampling\n');
            end
            %Sampling
            [stratifiedRankedSetSampleData,stratifiedRankedSetSampleVar(i,j),stratifiedRankedSetSampleStd(i,j),stratifiedRankedSetSampleMean(i,j)] = StratifiedRankedSetSampling(dataWithClassName, samplingLoop, stratifyNumber, sampleSize);
            
            %Clustering
            [stratifiedRankedSetSampleClusteringIndex,stratifiedRankedSetSampleClusteringCenter,stratifiedRankedSetSampleClusteringSumDistanceTemp] = clusteringFunction(stratifiedRankedSetSampleData(:,1),numberOfCluster);
            stratifiedRankedSetSampleClusteringSumDistance(i,j) = sum(stratifiedRankedSetSampleClusteringSumDistanceTemp);
            
            %Run Clustering Evaluation
            stratifiedRankedSetSampleClusterEval(i,j) = ClusteringEval(stratifiedRankedSetSampleData(:,2),stratifiedRankedSetSampleClusteringIndex,numberOfCluster);
            
        end
        
        if showMessage
            %fprintf(sprintf('Secend Loop %d End \n \n',j));
        end
    end
    
    
    
    if showMessage
        fprintf(sprintf('Main Loop %d End \n \n',i));
    end
end

%--------------------------------------------------------------------------
% Show Result
%--------------------------------------------------------------------------
if showMessage
    if runSimpleRandomSampling
        fprintf('Simple Random Sampling Summery:\n');
        %fprintf(sprintf('Avrage of Var %d \n',mean(simpleRandomSampleVar)));
        %fprintf(sprintf('Avrage of Std %d \n',mean(simpleRandomSampleStd)));
        %fprintf(sprintf('Avrage of Mean %d \n',mean(simpleRandomSampleMean)));
        %fprintf(sprintf('Avrage of Clustering Sum Distance %d \n',mean(simpleRandomSampleClusteringSumDistance)));
        fprintf(sprintf('Avrage of Clustering Eval %d \n',mean(reshape(simpleRandomSampleClusterEval,[],1))));
    end
    
    if runRankedSetSampling
        fprintf('Ranked Set Sampling Summery:\n');
        %fprintf(sprintf('Avrage of Var %d \n',mean(rankedSetSampleVar)));
        %fprintf(sprintf('Avrage of Std %d \n',mean(rankedSetSampleStd)));
        %fprintf(sprintf('Avrage of Mean %d \n',mean(rankedSetSampleMean)));
        %fprintf(sprintf('Avrage of Clustering Sum Distance %d \n',mean(rankedSetSampleClusteringSumDistance)));
        fprintf(sprintf('Avrage of Clustering Eval %d \n',mean(reshape(rankedSetSampleClusterEval,[],1))));
        
    end
    
    if runStratifiedRandomSampling
        fprintf('Stratified Random Sampling Summery:\n');
        %fprintf(sprintf('Avrage of Var %d \n',mean(stratifiedRandomSampleVar)));
        %fprintf(sprintf('Avrage of Std %d \n',mean(stratifiedRandomSampleStd)));
        %fprintf(sprintf('Avrage of Mean %d \n',mean(stratifiedRandomSampleMean)));
        %fprintf(sprintf('Avrage of Clustering Sum Distance %d \n',mean(stratifiedRandomSampleClusteringSumDistance)));
        fprintf(sprintf('Avrage of Clustering Eval %d \n',mean(reshape(stratifiedRandomSampleClusterEval,[],1))));
    end
    
    if runStratifiedRankedSetSampling
        fprintf('Stratified RankedSet Sampling Summery:\n');
        %fprintf(sprintf('Avrage of Var %d \n',mean(stratifiedRankedSetSampleVar)));
        %fprintf(sprintf('Avrage of Std %d \n',mean(stratifiedRankedSetSampleStd)));
        %fprintf(sprintf('Avrage of Mean %d \n',mean(stratifiedRankedSetSampleMean)));
        %fprintf(sprintf('Avrage of Clustering Sum Distance %d \n',mean(stratifiedRankedSetSampleClusteringSumDistance)));
        fprintf(sprintf('Avrage of Clustering Eval %d \n',mean(reshape(stratifiedRankedSetSampleClusterEval,[],1))));
        
    end
    if runAllData
        fprintf('All Data Summery:\n');
        %fprintf(sprintf('Avrage of Var %d \n',mean(stratifiedRankedSetSampleVar)));
        %fprintf(sprintf('Avrage of Std %d \n',mean(stratifiedRankedSetSampleStd)));
        %fprintf(sprintf('Avrage of Mean %d \n',mean(stratifiedRankedSetSampleMean)));
        %fprintf(sprintf('Avrage of Clustering Sum Distance %d \n',mean(stratifiedRankedSetSampleClusteringSumDistance)));
        fprintf(sprintf('Avrage of Clustering Eval %d \n',mean(allDataClusterEval)));
        
    end
end

% evalclusters(simpleRandomSampleData(:,1),clusteringFunction,'silhouette','klist',numberOfCluster)
% evalclusters(rankedSetSampleData(:,1),clusteringFunction,'silhouette','klist',numberOfCluster)
% evalclusters(stratifiedRandomSampleData(:,1),clusteringFunction,'silhouette','klist',numberOfCluster)
% evalclusters(stratifiedRankedSetSampleData(:,1),clusteringFunction,'silhouette','klist',numberOfCluster)
% evalclusters(dataWithClassName(:,1),clusteringFunction,'silhouette','klist',numberOfCluster)

%--------------------------------------------------------------------------
% End of Main Loop
%--------------------------------------------------------------------------
fprintf('End\n');