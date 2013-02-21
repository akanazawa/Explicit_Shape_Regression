function do_esr(config_file)
%%%%%%%%%%%%%%%%%%%%
% The driver code for the implementation of:
% Face Alignment by Explicit Shape Regression
% Xudong Cao, Yichen Wei, Fang Wen, and Jian Sun
% IEEE Conference on Computer Vision and Pattern Recognition
% (CVPR), 2012
%
% Angjoo Kanazawa July 12 2012
%%%%%%%%%%%%%%%%%%%%

%% load paramters
eval(config_file);
fprintf('starting ESR with N=%d A=%d T=%d K=%d P=%d beta=%d saving to %s\n',...
        param.N,param.A, param.T, param.K, param.P, param.beta, path.trainModel);
%% Get training data
if ~exist(path.trainData, 'file')
    [Is, Sgts] = getData(DIR, path.trainNames);
    save(path.trainData, 'Is', 'Sgts');
else
    fprintf('loading training data..\n');
    load(path.trainData);
end
% imageId = 1:size(Is,3);
imageId = 1:2000;
Is = Is(:,:,imageId);
Sgts = Sgts(:,:, imageId);

if param.N == 1000
    imageId = 1:2:2000;
    Is = Is(:,:, imageId);
    Sgts = Sgts(:,:, imageId);
end

%% Compute meanshape
[meanShape meanBox alignedGt] = computeMeanShape(Sgts, D, Is);

%% Train ferns
if ~exist(path.trainModel, 'file')
    %% Augment the training data 
    numTrain = size(Is, 3);
    Sgts = reshape(Sgts, 2*D.nParts, numTrain); % turn into vector
    alignedGt = reshape(alignedGt, 2*D.nParts, numTrain);
    initShapeAug = zeros(2*D.nParts, numTrain*param.A);
    SgtsAug = zeros(2*D.nParts, numTrain*param.A);
    augInds = zeros(1, numTrain*param.A);
    for i=1:numTrain
        % pick any shape other than this
        inds = randsample([1:i-1, i+1:numTrain], param.A);
        initShapeAug(:, (i-1)*param.A+1:i*param.A) = Sgts(:,inds);        
        SgtsAug(:, (i-1)*param.A+1:i*param.A) = repmat(Sgts(:, i), 1, param.A);
        augInds((i-1)*param.A+1:i*param.A) = i;
    end
    Sgts = [Sgts, SgtsAug];
    initShape = [repmat(meanShape(:), 1, numTrain), initShapeAug];
    [Rs, St]= do_train(Is, initShape, Sgts, meanShape, augInds, ...
                       meanBox, param, D); 
    save(path.trainModel, 'Rs', 'config_file', 'St');
    SgtsTrain = Sgts(:, 1:numTrain);
    clear Is Sgts;
else 
    load(path.trainModel); 
    fprintf('loading %s..\n', path.trainModel);
    numTrain = size(Sgts, 3);
    SgtsTrain = reshape(Sgts, 2*D.nParts, numTrain); % turn into vector
    alignedGt = reshape(alignedGt, 2*D.nParts, numTrain);
end

%% Get Test set
if ~exist(path.testData, 'file')
    [Is, Sgts] = getData(DIR, path.testNames);
    save(path.testData, 'Is', 'Sgts');
else
    fprintf('loading testing data..\n');
    load(path.testData);
    % fprintf('going with training data..\n');
    % numTrial = 1;
end


%% Run on Test
if param.N==2000
    numTest = 200;
else
    numTest = size(Is, 3);
end
imgInds = randperm(numTest);
Ndo = numTest;
StAll = zeros(2*D.nParts, Ndo);
initShapes = zeros(2*D.nParts, param.numTrials, Ndo);
score = zeros(Ndo,1);
scores = zeros(D.nParts,Ndo);
seeds = randi(numTrain, param.numTrials-1, Ndo);
VERBOSE = 0; 
fprintf('start testing on %d images\n', Ndo);
total = 0;
totalTime = 0;
for i = 1:Ndo
    ind = imgInds(i);
    seed = seeds(:, i);
    initShapes(:, :,i) = [meanShape(:), alignedGt(:, seed)];    
    % initShapes = [meanShape(:), SgtsTrain(:, randi(numTrain, param.numTrials-1, 1))];    
    [StAll(:,i), score(i)] = do_test(Rs, Is(:,:,ind), initShapes(:,:,i), ...
                                     Sgts(:,:,ind), meanShape, param ,D,VERBOSE); 
    
    interocular = sqrt(sum((Sgts(:, D.leye,ind) - Sgts(:, D.reye,ind)).^2));
    scores(:,i) = bsxfun(@rdivide, sqrt(sum((Sgts(:,:,ind) - reshape(StAll(:, i), [2 D.nParts])).^2)), interocular)';               
    VERBOSE = 0;
end
fprintf('average test error after %g trials %g, median:%g\n', param.numTrials,mean(scores(:)),median(scores(:)));
eval(config_file);
save(path.scores, 'score', 'StAll');

keyboard

[test, worstToBest] = sort(score, 'descend');

saveResults(Is(:,:,imgInds(worstToBest)), StAll(:, worstToBest), Sgts(:,:,imgInds(worstToBest)), initShapes(:,:,worstToBest),DIR, 'testimg',D);

% plots
allRange = zeros(param.T*2, numTest); 
for i=1:200
    allRange(:, i) = featureRange(Rs, Is(:,:, i), meanShape(:), Sgts(:,:,i), ...
                 meanShape, D, param);
end
% visualize
hfig = featureSelected(Rs, meanShape, meanBox, param, D);
% hfig = featureSelected_pair(Rs, meanShape, meanBox, param, D);
% range = featureRangeAll(Rs, meanShape, param, D);
sfigure; plot(1:param.T*2, mean(allRange, 2), '.-');
xlabel('T (in 1/2 increments)'); 
ylabel('avg range of selected Features');
plotpca(Rs, meanShape, param, D);
title(sprintf('PCA of Gaussian-Kmeans-SIFT N=%d T=%d K=%d A=%d b=%d',...
              param.N, param.T, param.K, param.A, param.beta));

title(sprintf('PCA of Gaussian-Kmeans-KNN=3 N=%d T=%d K=%d A=%d b=%d',...
              param.N, param.T, param.K, param.A, param.beta));

