clc
clear all;
%--------------------------------------------------------------------------
% Setting
%-------------------------------------------------------------------------
initPopulation = 100000;% Number of Population
sampleSizeRenge = [5 10 20];% Size of Sample
samplingLoopRenge= [20 10 5];% Number of Loop for Sampling
% stratifyNumber = samplingLoop;
firstMeanStart = -4;% Mean
meanIntervalRenge = [1 2 3 4];% Interval between of 2 Mean
SIGMARenge =[1 2 3];% Sigma
numberOfCluster = 5; % Number Of cluster for init(generate data) and k-means
mainLoop = 10; % Number of Loop for All Program
secendLoop = 20; % Number of Loop for each main loop
% Clustering Algorithm
clusteringFunction = @(X,K) kmeans(X, K, ...
    'EmptyAction','singleton', 'Replicates',5, ...
    'Distance','sqeuclidean','MaxIter',500,...
    'Display','off','Options',statset('UseParallel',0));
showMessage = 1; % On/Off Command Window Message
showSubMessage = 0;% On/Off Command Window Sub Message
showVarMessage = 0;% On/Off Command Window Var Message
showStdMessage = 0;% On/Off Command Window Std Message
showMeanMessage = 0;% On/Off Command Window Mean Message
showClusterDistMessage = 0;% On/Off Command Window ClusterDist Message
showClusterEvalMessage = 1;% On/Off Command Window ClusterEval Message

% On/Off Run Clustering on each Sampling Method:
runAllData = 0; % Run Clustering on All Data

runSimpleRandomSampling = 1; % Run Simple Random

runRankedSetSamplingDiagonal = 1; % Run Ranked Set(type:Diagonal)
runRankedSetSamplingMiddel = 0; % Run Ranked Set(type:Middel)
runRankedSetSamplingFirst = 0; % Run Ranked Set(type:First)
runRankedSetSamplingLast = 0; % Run Ranked Set(type:Last)

runStratifiedRandomSampling = 1; % Run Stratified Random

runStratifiedRankedSetSamplingDiagonal = 1; % Run Stratified Ranked Set (type:Diagonal)
runStratifiedRankedSetSamplingMiddel = 1; % Run Stratified Ranked Set (type:Middel)
runStratifiedRankedSetSamplingFirst = 0; % Run Stratified Ranked Set (type:First)
runStratifiedRankedSetSamplingLast = 0; % Run Stratified Ranked Set (type:Last)

%--------------------------------------------------------------------------
% Start Main
%--------------------------------------------------------------------------
fprintf('Start\n');

%--------------------------------------------------------------------------
% Create Zero Vector For SpeedUp
%--------------------------------------------------------------------------
% simpleRandomSampleVar = zeros(mainLoop,secendLoop);
% simpleRandomSampleStd = zeros(mainLoop,secendLoop);
% simpleRandomSampleMean = zeros(mainLoop,secendLoop);
% simpleRandomSampleClusteringSumDistance = zeros(mainLoop,secendLoop);
% simpleRandomSampleClusterEval = zeros(mainLoop,secendLoop);

% rankedSetSampleVar = zeros(mainLoop,secendLoop);
% rankedSetSampleStd = zeros(mainLoop,secendLoop);
% rankedSetSampleMean = zeros(mainLoop,secendLoop);
% rankedSetSampleClusteringSumDistance = zeros(mainLoop,secendLoop);
% rankedSetSampleClusterEval = zeros(mainLoop,secendLoop);

% stratifiedRandomSampleVar = zeros(mainLoop,secendLoop);
% stratifiedRandomSampleStd = zeros(mainLoop,secendLoop);
% stratifiedRandomSampleMean = zeros(mainLoop,secendLoop);
% stratifiedRandomSampleClusteringSumDistance = zeros(mainLoop,secendLoop);
% stratifiedRandomSampleClusterEval = zeros(mainLoop,secendLoop);

% stratifiedRankedSetSampleVar = zeros(mainLoop,secendLoop);
% stratifiedRankedSetSampleStd = zeros(mainLoop,secendLoop);
% stratifiedRankedSetSampleMean = zeros(mainLoop,secendLoop);
% stratifiedRankedSetSampleClusteringSumDistance = zeros(mainLoop,secendLoop);
% stratifiedRankedSetSampleClusterEval = zeros(mainLoop,secendLoop);
if runAllData
    allDataClusteringSumDistance = zeros(mainLoop,1);
    allDataClusterEval = zeros(mainLoop,1);
end
% Why secendLoop * 10 ?
% 10 for SimpleRandomSampling + RankedSetSamplingDiagonal +
%        RankedSetSamplingMiddel + RankedSetSamplingFirst +
%        RankedSetSamplingLast + StratifiedRandomSampling +
%        StratifiedRankedSetSamplingDiagonal + StratifiedRankedSetSamplingMiddel +
%        StratifiedRankedSetSamplingFirst + StratifiedRankedSetSamplingLast

samplingVar = zeros(mainLoop,secendLoop*10);
samplingStd = zeros(mainLoop,secendLoop*10);
samplingMean = zeros(mainLoop,secendLoop*10);
samplingClusterDist = zeros(mainLoop,secendLoop*10);
samplingClusterEval = zeros(mainLoop,secendLoop*10);

resultRowCount = length(sampleSizeRenge)*length(meanIntervalRenge)*length(SIGMARenge);
resultClusterEval = zeros(resultRowCount,15);

%--------------------------------------------------------------------------
appLoop = 0;
for iSampleSize = 1 : length(sampleSizeRenge)
        for iMeanInterval = 1 : length(meanIntervalRenge)
            for iSIGMA = 1 : length(SIGMARenge)
                
                appLoop = appLoop +1; 
                
                sampleSize = sampleSizeRenge(iSampleSize);
                samplingLoop = samplingLoopRenge(iSampleSize);
                stratifyNumber = samplingLoop;
                meanInterval = meanIntervalRenge(iMeanInterval);
                SIGMA = SIGMARenge(iSIGMA);
                if showMessage
                    fprintf(sprintf('\n --------- %d ----------- \n',appLoop));
                    fprintf(sprintf('sampleSize = %d \n',sampleSize));
                    fprintf(sprintf('samplingLoop = %d \n',samplingLoop));
                    fprintf(sprintf('meanInterval = %d \n',meanInterval));
                    fprintf(sprintf('SIGMA = %d \n',SIGMA));
                end
                %--------------------------------------------------------------------------
                %--------------------------------------------------------------------------
                %--------------------------------------------------------------------------
                %--------------------------------------------------------------------------
                
                
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
                    %--------------------------------------------------------------------------
                    % Start of Main Loop
                    %--------------------------------------------------------------------------
                    %--------------------------------------------------------------------------
                    for j = 1: secendLoop
                        if showSubMessage
                            fprintf(sprintf('Secend Loop %d Start \n',j));
                        end
                        
                        %----------------------------------------------------------------------
                        %----------------------------------------------------------------------
                        %Simple Random Sampling
                        %----------------------------------------------------------------------
                        %----------------------------------------------------------------------
                        if runSimpleRandomSampling
                            if showSubMessage
                                fprintf('Simple Random Sampling\n');
                            end
                            %Sampling
                            [simpleRandomSampleData,samplingVar(i,j),samplingStd(i,j),samplingMean(i,j)] = SimpleRandomSampling(dataWithClassName, samplingLoop, initPopulation, sampleSize);
                            
                            %Clustering
                            [simpleRandomSampleClusteringIndex,simpleRandomSampleClusteringCenter,simpleRandomSampleClusteringSumDistanceTemp] = clusteringFunction(simpleRandomSampleData(:,1),numberOfCluster);
                            samplingClusterDist(i,j) = sum(simpleRandomSampleClusteringSumDistanceTemp);
                            
                            %Clustering Evaluation
                            samplingClusterEval(i,j) = ClusteringEval(simpleRandomSampleData(:,2),simpleRandomSampleClusteringIndex,numberOfCluster);
                        end
                        
                        %----------------------------------------------------------------------
                        %----------------------------------------------------------------------
                        %Ranked Set Sampling
                        %----------------------------------------------------------------------
                        %----------------------------------------------------------------------
                        
                        %----------------------------------------------------------------------
                        %Diagonal
                        %----------------------------------------------------------------------
                        if runRankedSetSamplingDiagonal
                            if showSubMessage
                                fprintf('Ranked Set Sampling\n');
                            end
                            %Sampling
                            [rankedSetSampleData,samplingVar(i,j+secendLoop),samplingStd(i,j+secendLoop),samplingMean(i,j+secendLoop)] = RankedSetSampling(dataWithClassName, samplingLoop, initPopulation, sampleSize,1);
                            
                            %Clustering
                            [rankedSetSampleClusteringIndex,rankedSetSampleClusteringCenter,rankedSetSampleClusteringSumDistanceTemp] = clusteringFunction(rankedSetSampleData(:,1),numberOfCluster);
                            samplingClusterDist(i,j+secendLoop) = sum(rankedSetSampleClusteringSumDistanceTemp);
                            
                            %Clustering Evaluation
                            samplingClusterEval(i,j+secendLoop) = ClusteringEval(rankedSetSampleData(:,2),rankedSetSampleClusteringIndex,numberOfCluster);
                            
                        end
                        
                        %----------------------------------------------------------------------
                        %Middel
                        %----------------------------------------------------------------------
                        if runRankedSetSamplingMiddel
                            if showSubMessage
                                fprintf('Ranked Set Sampling\n');
                            end
                            %Sampling
                            [rankedSetSampleData,samplingVar(i,j+(secendLoop*2)),samplingStd(i,j+(secendLoop*2)),samplingMean(i,j+(secendLoop*2))] = RankedSetSampling(dataWithClassName, samplingLoop, initPopulation, sampleSize,2);
                            
                            %Clustering
                            [rankedSetSampleClusteringIndex,rankedSetSampleClusteringCenter,rankedSetSampleClusteringSumDistanceTemp] = clusteringFunction(rankedSetSampleData(:,1),numberOfCluster);
                            samplingClusterDist(i,j+(secendLoop*2)) = sum(rankedSetSampleClusteringSumDistanceTemp);
                            
                            %Clustering Evaluation
                            samplingClusterEval(i,j+(secendLoop*2)) = ClusteringEval(rankedSetSampleData(:,2),rankedSetSampleClusteringIndex,numberOfCluster);
                            
                        end
                        
                        %----------------------------------------------------------------------
                        %First Point
                        %----------------------------------------------------------------------
                        if runRankedSetSamplingFirst
                            if showSubMessage
                                fprintf('Ranked Set Sampling\n');
                            end
                            %Sampling
                            [rankedSetSampleData,samplingVar(i,j+(secendLoop*3)),samplingStd(i,j+(secendLoop*3)),samplingMean(i,j+(secendLoop*3))] = RankedSetSampling(dataWithClassName, samplingLoop, initPopulation, sampleSize,3);
                            
                            %Clustering
                            [rankedSetSampleClusteringIndex,rankedSetSampleClusteringCenter,rankedSetSampleClusteringSumDistanceTemp] = clusteringFunction(rankedSetSampleData(:,1),numberOfCluster);
                            samplingClusterDist(i,j+(secendLoop*3)) = sum(rankedSetSampleClusteringSumDistanceTemp);
                            
                            %Clustering Evaluation
                            samplingClusterEval(i,j+(secendLoop*3)) = ClusteringEval(rankedSetSampleData(:,2),rankedSetSampleClusteringIndex,numberOfCluster);
                            
                        end
                        
                        %----------------------------------------------------------------------
                        %Last Point
                        %----------------------------------------------------------------------
                        if runRankedSetSamplingLast
                            if showSubMessage
                                fprintf('Ranked Set Sampling\n');
                            end
                            %Sampling
                            [rankedSetSampleData,samplingVar(i,j+(secendLoop*4)),samplingStd(i,j+(secendLoop*4)),samplingMean(i,j+(secendLoop*4))] = RankedSetSampling(dataWithClassName, samplingLoop, initPopulation, sampleSize,4);
                            
                            %Clustering
                            [rankedSetSampleClusteringIndex,rankedSetSampleClusteringCenter,rankedSetSampleClusteringSumDistanceTemp] = clusteringFunction(rankedSetSampleData(:,1),numberOfCluster);
                            samplingClusterDist(i,j+(secendLoop*4)) = sum(rankedSetSampleClusteringSumDistanceTemp);
                            
                            %Clustering Evaluation
                            samplingClusterEval(i,j+(secendLoop*4)) = ClusteringEval(rankedSetSampleData(:,2),rankedSetSampleClusteringIndex,numberOfCluster);
                            
                        end
                        
                        %----------------------------------------------------------------------
                        %----------------------------------------------------------------------
                        %Stratified Random Sampling
                        %----------------------------------------------------------------------
                        %----------------------------------------------------------------------
                        if runStratifiedRandomSampling
                            if showSubMessage
                                fprintf('Stratified Random Sampling\n');
                            end
                            %Sampling
                            [stratifiedRandomSampleData,samplingVar(i,j+(secendLoop*5)),samplingStd(i,j+(secendLoop*5)),samplingMean(i,j+(secendLoop*5))] = StratifiedRandomSampling(dataWithClassName, samplingLoop, stratifyNumber, sampleSize);
                            
                            %Clustering
                            [stratifiedRandomSampleClusteringIndex,stratifiedRandomSampleClusteringCenter,stratifiedRandomSampleClusteringSumDistanceTemp] = clusteringFunction(stratifiedRandomSampleData(:,1),numberOfCluster);
                            samplingClusterDist(i,j+(secendLoop*5)) = sum(stratifiedRandomSampleClusteringSumDistanceTemp);
                            
                            %Run Clustering Evaluation
                            samplingClusterEval(i,j+(secendLoop*5)) = ClusteringEval(stratifiedRandomSampleData(:,2),stratifiedRandomSampleClusteringIndex,numberOfCluster);
                            
                        end
                        
                        %----------------------------------------------------------------------
                        %----------------------------------------------------------------------
                        %Stratified Ranked Set Sampling
                        %----------------------------------------------------------------------
                        %----------------------------------------------------------------------
                        
                        %----------------------------------------------------------------------
                        %Diagonal
                        %----------------------------------------------------------------------
                        if runStratifiedRankedSetSamplingDiagonal
                            if showSubMessage
                                fprintf('Stratified RankedSet Sampling\n');
                            end
                            %Sampling
                            [stratifiedRankedSetSampleData,samplingVar(i,j+(secendLoop*6)),samplingStd(i,j+(secendLoop*6)),samplingMean(i,j+(secendLoop*6))] = StratifiedRankedSetSampling(dataWithClassName, samplingLoop, stratifyNumber, sampleSize,1);
                            
                            %Clustering
                            [stratifiedRankedSetSampleClusteringIndex,stratifiedRankedSetSampleClusteringCenter,stratifiedRankedSetSampleClusteringSumDistanceTemp] = clusteringFunction(stratifiedRankedSetSampleData(:,1),numberOfCluster);
                            samplingClusterDist(i,j+(secendLoop*6)) = sum(stratifiedRankedSetSampleClusteringSumDistanceTemp);
                            
                            %Run Clustering Evaluation
                            samplingClusterEval(i,j+(secendLoop*6)) = ClusteringEval(stratifiedRankedSetSampleData(:,2),stratifiedRankedSetSampleClusteringIndex,numberOfCluster);
                            
                        end
                        
                        %----------------------------------------------------------------------
                        %Middel
                        %----------------------------------------------------------------------
                        if runStratifiedRankedSetSamplingMiddel
                            if showSubMessage
                                fprintf('Stratified RankedSet Sampling\n');
                            end
                            %Sampling
                            [stratifiedRankedSetSampleData,samplingVar(i,j+(secendLoop*7)),samplingStd(i,j+(secendLoop*7)),samplingMean(i,j+(secendLoop*7))] = StratifiedRankedSetSampling(dataWithClassName, samplingLoop, stratifyNumber, sampleSize,2);
                            
                            %Clustering
                            [stratifiedRankedSetSampleClusteringIndex,stratifiedRankedSetSampleClusteringCenter,stratifiedRankedSetSampleClusteringSumDistanceTemp] = clusteringFunction(stratifiedRankedSetSampleData(:,1),numberOfCluster);
                            samplingClusterDist(i,j+(secendLoop*7)) = sum(stratifiedRankedSetSampleClusteringSumDistanceTemp);
                            
                            %Run Clustering Evaluation
                            samplingClusterEval(i,j+(secendLoop*7)) = ClusteringEval(stratifiedRankedSetSampleData(:,2),stratifiedRankedSetSampleClusteringIndex,numberOfCluster);
                            
                        end
                        
                        %----------------------------------------------------------------------
                        %First Point
                        %----------------------------------------------------------------------
                        if runStratifiedRankedSetSamplingFirst
                            if showSubMessage
                                fprintf('Stratified RankedSet Sampling\n');
                            end
                            %Sampling
                            [stratifiedRankedSetSampleData,samplingVar(i,j+(secendLoop*8)),samplingStd(i,j+(secendLoop*8)),samplingMean(i,j+(secendLoop*8))] = StratifiedRankedSetSampling(dataWithClassName, samplingLoop, stratifyNumber, sampleSize,3);
                            
                            %Clustering
                            [stratifiedRankedSetSampleClusteringIndex,stratifiedRankedSetSampleClusteringCenter,stratifiedRankedSetSampleClusteringSumDistanceTemp] = clusteringFunction(stratifiedRankedSetSampleData(:,1),numberOfCluster);
                            samplingClusterDist(i,j+(secendLoop*8)) = sum(stratifiedRankedSetSampleClusteringSumDistanceTemp);
                            
                            %Run Clustering Evaluation
                            samplingClusterEval(i,j+(secendLoop*8)) = ClusteringEval(stratifiedRankedSetSampleData(:,2),stratifiedRankedSetSampleClusteringIndex,numberOfCluster);
                            
                        end
                        
                        %----------------------------------------------------------------------
                        %Last Point
                        %----------------------------------------------------------------------
                        if runStratifiedRankedSetSamplingLast
                            if showSubMessage
                                fprintf('Stratified RankedSet Sampling\n');
                            end
                            %Sampling
                            [stratifiedRankedSetSampleData,samplingVar(i,j+(secendLoop*9)),samplingStd(i,j+(secendLoop*9)),samplingMean(i,j+(secendLoop*9))] = StratifiedRankedSetSampling(dataWithClassName, samplingLoop, stratifyNumber, sampleSize,4);
                            
                            %Clustering
                            [stratifiedRankedSetSampleClusteringIndex,stratifiedRankedSetSampleClusteringCenter,stratifiedRankedSetSampleClusteringSumDistanceTemp] = clusteringFunction(stratifiedRankedSetSampleData(:,1),numberOfCluster);
                            samplingClusterDist(i,j+(secendLoop*9)) = sum(stratifiedRankedSetSampleClusteringSumDistanceTemp);
                            
                            %Run Clustering Evaluation
                            samplingClusterEval(i,j+(secendLoop*9)) = ClusteringEval(stratifiedRankedSetSampleData(:,2),stratifiedRankedSetSampleClusteringIndex,numberOfCluster);
                            
                        end
                        
                        if showSubMessage
                            fprintf(sprintf('Secend Loop %d End \n \n',j));
                        end
                    end
                    
                    
                    
                    if showSubMessage
                        fprintf(sprintf('Main Loop %d End \n \n',i));
                    end
                end
                
                %--------------------------------------------------------------------------
                % Show Result
                %--------------------------------------------------------------------------
                
                

                if showVarMessage
                    % TODO!
                    fprintf('Simple Random Sampling Var:');
                    fprintf(sprintf('Avrage of Clustering Eval %d \n',mean(reshape(simpleRandomSampleClusterEval,[],1))));
                end
                
                if showStdMessage
                    % TODO!
                    fprintf('Simple Random Sampling Var:');
                    fprintf(sprintf('Avrage of Clustering Eval %d \n',mean(reshape(simpleRandomSampleClusterEval,[],1))));
                end
                
                if showMeanMessage
                    % TODO!
                    fprintf('Simple Random Sampling Var:');
                    fprintf(sprintf('Avrage of Clustering Eval %d \n',mean(reshape(simpleRandomSampleClusterEval,[],1))));
                end
                
                if showClusterDistMessage
                    % TODO!
                    fprintf('Simple Random Sampling Var:');
                    fprintf(sprintf('Avrage of Clustering Eval %d \n',mean(reshape(simpleRandomSampleClusterEval,[],1))));
                end
                
                resultClusterEval(appLoop,1) = sampleSize;
                resultClusterEval(appLoop,2) = samplingLoop;
                resultClusterEval(appLoop,3) = meanInterval;
                resultClusterEval(appLoop,4) = SIGMA;
                
                SRStmp = mean(reshape(samplingClusterEval(:,1:secendLoop),[],1));
                resultClusterEval(appLoop,5) = SRStmp;
                
                RSSDtmp = mean(reshape(samplingClusterEval(:,secendLoop+1:secendLoop*2),[],1));
                resultClusterEval(appLoop,6) = RSSDtmp;
                
                RSSMtmp = mean(reshape(samplingClusterEval(:,(secendLoop*2)+1:secendLoop*3),[],1));
                resultClusterEval(appLoop,7) = RSSMtmp;
                
                RSSFtmp = mean(reshape(samplingClusterEval(:,(secendLoop*3)+1:secendLoop*4),[],1));
                resultClusterEval(appLoop,8) = RSSFtmp;
                
                RSSLtmp = mean(reshape(samplingClusterEval(:,(secendLoop*4)+1:secendLoop*5),[],1));
                resultClusterEval(appLoop,9) = RSSLtmp;
                
                StrafRStmp = mean(reshape(samplingClusterEval(:,(secendLoop*5)+1:secendLoop*6),[],1));
                resultClusterEval(appLoop,10) = StrafRStmp;
                
                StrafRSSDtmp = mean(reshape(samplingClusterEval(:,(secendLoop*6)+1:secendLoop*7),[],1));
                resultClusterEval(appLoop,11) = StrafRSSDtmp;
                
                StrafRSSMtmp = mean(reshape(samplingClusterEval(:,(secendLoop*7)+1:secendLoop*8),[],1));
                resultClusterEval(appLoop,12) = StrafRSSMtmp;
                
                StrafRSSFtmp = mean(reshape(samplingClusterEval(:,(secendLoop*8)+1:secendLoop*9),[],1));
                resultClusterEval(appLoop,13) = StrafRSSFtmp;
                
                StrafRSSLtmp = mean(reshape(samplingClusterEval(:,(secendLoop*9)+1:secendLoop*10),[],1));
                resultClusterEval(appLoop,14) = StrafRSSLtmp;
                
                
                xlswrite('resultClusterEval.xlsx',resultClusterEval,1);
                
                if showClusterEvalMessage
                    fprintf(sprintf('SRS Clustering Eval: %d \n',SRStmp));
                    fprintf(sprintf('RSS Diagonal Clustering Eval: %d \n',RSSDtmp));
                    fprintf(sprintf('RSS Middel Clustering Eval: %d \n',RSSMtmp));
                    fprintf(sprintf('RSS First Clustering Eval: %d \n',RSSFtmp));
                    fprintf(sprintf('RSS Last Clustering Eval: %d \n',RSSLtmp));
                    fprintf(sprintf('Strafied RS Clustering Eval: %d \n',StrafRStmp));
                    fprintf(sprintf('Strafied RRS Diagonal Clustering Eval: %d \n',StrafRSSDtmp));
                    fprintf(sprintf('Strafied RRS Middel Clustering Eval: %d \n',StrafRSSMtmp));
                    fprintf(sprintf('Strafied RRS First Clustering Eval: %d \n',StrafRSSFtmp));
                    fprintf(sprintf('Strafied RRS Last Clustering Eval: %d \n',StrafRSSLtmp));
                end
                
                if runAllData
                    fprintf('All Data Summery:\n');
                    %fprintf(sprintf('Avrage of Var %d \n',mean(stratifiedRankedSetSampleVar)));
                    %fprintf(sprintf('Avrage of Std %d \n',mean(stratifiedRankedSetSampleStd)));
                    %fprintf(sprintf('Avrage of Mean %d \n',mean(stratifiedRankedSetSampleMean)));
                    %fprintf(sprintf('Avrage of Clustering Sum Distance %d \n',mean(stratifiedRankedSetSampleClusteringSumDistance)));
                    fprintf(sprintf('Avrage of Clustering Eval %d \n',mean(allDataClusterEval)));
                    
                end
                
                
                %--------------------------------------------------------------------------
                %--------------------------------------------------------------------------
                %--------------------------------------------------------------------------
                %--------------------------------------------------------------------------
                
                
            end
        end
    
end


%--------------------------------------------------------------------------


% evalclusters(simpleRandomSampleData(:,1),clusteringFunction,'silhouette','klist',numberOfCluster)
% evalclusters(rankedSetSampleData(:,1),clusteringFunction,'silhouette','klist',numberOfCluster)
% evalclusters(stratifiedRandomSampleData(:,1),clusteringFunction,'silhouette','klist',numberOfCluster)
% evalclusters(stratifiedRankedSetSampleData(:,1),clusteringFunction,'silhouette','klist',numberOfCluster)
% evalclusters(dataWithClassName(:,1),clusteringFunction,'silhouette','klist',numberOfCluster)

%--------------------------------------------------------------------------
% End of Main Loop
%--------------------------------------------------------------------------
fprintf('End\n');