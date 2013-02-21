%%%%%%%%%%%%%%%%%%%%
% Configuration file for the implementation of Face Alignment by
% Explicit Shape Regression
% to allow easier replication of experiments with various
% parameters
%%%%%%%%%%%%%%%%%%%%
s = RandStream('mcg16807','Seed',0);
RandStream.setGlobalStream(s);    
addpath('misc');

%% PATHS
DIR.img = '';
DIR.result = 'results/msra';
DIR.data = 'data/'; % save intermediate files

%path to the data file
path.trainData = [DIR.data, 'trainData.mat'];
% path.trainNames = 'peopleDevTrain.txt';
path.testData = [DIR.data, 'testData.mat'];
% path.testNames = '.peopleDevTest.txt';

%% ALGORITHM PARAMETERS
param.N = 2000; % number of images
param.P = 400;
param.T = 10;
param.K = 500; 
param.F = 5; 
param.numTrials = 5;
param.mode = 'original';

param.A = 20; 
param.beta = 1000;(param.N*param.A)/10;
path.trainModel = [DIR.data, sprintf('trainedFerns_default_B%d.mat', ...
                                     param.beta)];
% end

%% Dataset information
D.leye = 0; %left pupil ID
D.reye = 0; %right pupil ID
D.nParts = 0;
D.nRow = 250; D.nCol = 250; % all images are 250x250
% for plotting
D.connectedParts = {
        [1,2,4,6,8,7,5,3],
        [9 10 11],
        [14 13]};
