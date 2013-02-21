function [Rs, St] = do_train(Is, initShape, Sgts, meanShape, augInds, ...
                             meanBox, param, D)
%%%%%%%%%%%%%%%%%%%%
% Training of the ferns
%%%%%%%%%%%%%%%%%%%%
numImage = size(Is, 3);
numTrain = size(initShape, 2);
Is = reshape(Is, [D.nRow*D.nCol, numImage]);
offset = 0:(D.nRow*D.nCol):(numImage-1)*D.nRow*D.nCol;
offsetAug = (augInds-1).*D.nRow*D.nCol; 
St = initShape;
SgtsSub = reshape(Sgts, [2 D.nParts numTrain]);
allInterocular = sqrt(sum((SgtsSub(:, D.leye, :) - SgtsSub(:, D.reye, :)).^2));
[X Y] = meshgrid(1:param.P);
Rs(param.T) = struct('ferns', struct([]));
for t=1:param.T
    StSub = reshape(St, [2 D.nParts numTrain]);
    ME = bsxfun(@rdivide, sqrt(sum((SgtsSub - StSub).^2)), allInterocular);
    fprintf('total error @ t=%d:%g\n', t, mean(ME(:)));
    % sample P points wout replacement in the meanshape coordinate
    % sfigure; rectangle('Position', meanBox, 'EdgeColor', 'g', 'LineWidth', 3);    
    % hold on; axis image ij;
    % plot(xi, yi, 'c.') 
    %%index in (x,y) because shapes are in (x,y), not (r,c)=(y,x)
    % pinds = randsample(meanBox(4)*meanBox(3), param.P,1); %% OLD WAY
    pinds = randsample(meanBox(4)*meanBox(3), param.P); %% NEW WAY
    yi = rem(pinds-1, meanBox(4))+1 + meanBox(2); 
    xi = (pinds-yi)/meanBox(4) +1 + meanBox(1);

    % find the nearest landmark for each sample point
    dists = sqrt(sum(bsxfun(@minus, permute([xi, yi]', [1 3 2]), meanShape).^2, 1));
    [~, nearestInd] = min(permute(dists, [2 3 1]), [], 1);
    % compute the local coordinate wrt to the nearest landmark in meanshape
    localPts = [xi, yi]' - meanShape(:, nearestInd);

    %% transform current shape to mean shape, sample index
    localPtsTransform = zeros([size(localPts) numTrain]);
    for i=1:numTrain
        T = bestfit_nonreflective_similarity(StSub(:,:,i)',meanShape');
        localPtsTransform(:,:,i) = T'*[localPts; zeros(1, param.P)];
    end
    allSamplePts = round(StSub(:, nearestInd, :) + localPtsTransform);
    allSamplePts(1,allSamplePts(1,:,:) > D.nCol) = D.nCol;
    allSamplePts(2,allSamplePts(2,:,:) > D.nRow) = D.nRow;
    allSamplePts(allSamplePts < 1) = 1;
    % allInds = sub2ind([D.nRow, D.nCol], allSamplePts(2,:,:), allSamplePts(1,:,:));
    allInds = allSamplePts(2,:,:) + (allSamplePts(1,:,:)-1).*D.nRow;

    allInds = permute(allInds, [2 3 1]);
    allIndsOrig = bsxfun(@plus, allInds(:, 1:numImage), offset);
    allIndsAug =  bsxfun(@plus, allInds(:, numImage+1:end), offsetAug);
    allPs  = double([Is(allIndsOrig), Is(allIndsAug)]);
    
    % vis = randi(numImage, 1,1); 68610;%
    % visualize(reshape(Is(:, vis), [250 250]), reshape(initShape(:,vis), [2 55]),SgtsSub(:,:,vis), D, allSamplePts(:,:,vis), allSamplePts(:,:,vis))

    %% learn K many ferns
    allFern(param.K) = struct('fern', struct([]), 'output', []);
    % compute std(fi-fj) = sqrt((var(fi)+var(fj)-2*covFij(i,j)+eps));

    covFij = cov(allPs'); 
    varFij = diag(covFij);
    stdFij = reshape(varFij(X,:) + varFij(Y,:), param.P, param.P)-2.*covFij;
    stdFij = sqrt(stdFij+eps);
    StSubOld = StSub;
    nomeanAllPs = bsxfun(@minus,allPs,sum(allPs, 2)./numTrain);    
    for k=1:param.K
        deltaShapes = Sgts - St;
        fern(param.F) = struct('fi', [], 'fj', [], 'tau', 0);
        signature = zeros(param.F, numTrain);
        for f = 1:param.F        % learn F features
            % pick a random direction: sample from a zero mean gaussian
            direction = randn(2*D.nParts, 1);
            direction = direction./sqrt(sum(direction.^2));   % normalize it to unit vector
            y = sum(bsxfun(@times, deltaShapes, direction)); %project
            % sample covariance between Y and P:
            % sum((y-mean(y)).*(allPs(i,:)-mean(allPs(i,:))))/(numTrain-1);
            % correlation coefficient: corr(Y,Fi-Fj) = cov(Y, Fi-Fj)./std(Y)std(Fi-Fj)  
            covYF = sum(bsxfun(@times,(y-mean(y)),nomeanAllPs), 2)./(numTrain-1); 
            corr = reshape(covYF(Y,:) - covYF(X,:), param.P, param.P);
            corr = corr./(std(y).*(stdFij));             
            [maxVal, maxInd] = max(corr(:));  % find the most correlated
            [fi fj] = ind2sub(size(corr), maxInd);
            assert(corr(fi,fj) - maxVal < 1e-10)
            % to do randmoly
            % fi = randsample(param.P, 1); fj= randsample([1:fi, fi+1:param.P], 1);
            vals = allPs(fi,:) - allPs(fj,:);
            % maxVal  = max(vals).*.8; minVal = min(vals).*.8;
            % thresh = (maxVal-minVal)*rand(1)+minVal;
            thresh = randn(1).*std(vals)+mean(vals);

            % build ferns
            fern(f).fi = [localPts(:, fi); nearestInd(fi)];
            fern(f).fj = [localPts(:, fj); nearestInd(fj)];
            fern(f).tau = thresh;
            signature(f, vals > thresh) = 1;% mark the training data
        end
        %% divide into 2^F bins and update
        binAsg = sum(bsxfun(@times, signature, 2.^[0:param.F-1]'));
        deltaSb = zeros(2*D.nParts, 2.^param.F);        

        for b = 1:2^param.F
            nBin = sum(binAsg==b-1);
            if nBin ~= 0
                delta = sum(Sgts(:, binAsg==b-1) - St(:, binAsg==b-1), 2)./...
                        (nBin + param.beta);       
                St(:, binAsg==b-1) = bsxfun(@plus, St(:, binAsg==b-1), delta);
                deltaSb(:, b) = delta;
            end
        end
        allFern(k).fern = fern;
        allFern(k).output = deltaSb;        
        if mod(k, param.K/5)==0
            StSub = reshape(St, [2 D.nParts numTrain]);
            ME = bsxfun(@rdivide, sqrt(sum((SgtsSub - StSub).^2)), allInterocular);
            fprintf('trained %d/%d fern error:%g\n', k, param.K, mean(ME(:)));
        end
    end       %% end of K
    Rs(t).ferns = allFern;
end

