function [Is, Sgts] = getData(DIR, trainNames)
%%%%%%%%%%%%%%%%%%%%
% Extracting training images/ground truth parts 
%%%%%%%%%%%%%%%%%%%%

fid = fopen(trainNames);
numTrain = fscanf(fid, '%d', 1);
Is = uint8(zeros(D.nRow, D.nCol, numTrain));
Sgts = zeros(2, D.nParts, numTrain); % ignore the blocked for the moment
for i=1:numTrain
    name = fscanf(fid, '%s', 1);
    numImg = fscanf(fid, '%d', 1);
    fname = fullfile(DIR.img, name, sprintf('%s_%04d.jpg', name,numImg));
    fprintf('%d/%d @ %s..\n', i, numTrain, fname);
    matname = fullfile(DIR.img,name, sprintf('%s_%04d.mat', name, numImg));
    I = rgb2gray(imread(fname));
    % Is(:, :, :, i) = I;
    Is(:, :, i) = I;
    load(matname); % load parts
    Sgts(:, :, i) = parts(1:2, :);
    clear parts;
end
fclose(fid);

