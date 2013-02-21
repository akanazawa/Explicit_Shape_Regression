function [StMedian, fME] = do_test(Rs, I, initShapes, Sgt, meanShape, ...
                             param,D, VERBOSE)
%%%%%%%%%%%%%%%%%%%%
% Test on an single image
%%%%%%%%%%%%%%%%%%%%
interocular = sqrt(sum((Sgt(:, D.leye) - Sgt(:, D.reye)).^2));
offset = [0:2^param.F:(param.K-1)*2^param.F];
homoT = zeros(1,param.K*param.F);
StAll = zeros(2*D.nParts, param.numTrials);
for trial = 1:param.numTrials
    St = initShapes(:, trial);    
    % if VERBOSE,  sfigure; end
    for t=1:param.T
        ferns = [Rs(t).ferns.fern]; % (K*F)x1 struct
        allFi = [ferns.fi];
        allFj = [ferns.fj];    
        StSub = reshape(St, [2 D.nParts]);
        T = bestfit_nonreflective_similarity(StSub',meanShape');
        localPtsTransformFI = T'*[allFi(1:2,:); homoT];
        localPtsTransformFJ = T'*[allFj(1:2,:); homoT];    

        feat1 = round(StSub(:, allFi(3,:)) + localPtsTransformFI);
        feat1(1,feat1(1,:) > D.nCol) = D.nCol;
        feat1(2,feat1(2,:) > D.nRow) = D.nRow;
        feat1(feat1 < 1) = 1;
        feat2 = round(StSub(:, allFj(3,:)) + localPtsTransformFJ); 
        feat2(1,feat2(1,:) > D.nCol) = D.nCol;
        feat2(2,feat2(2,:) > D.nRow) = D.nRow;
        feat2(feat2 < 1) = 1;        

        allIndsFI = feat1(2,:,:) + (feat1(1,:,:)-1).*D.nRow;
        fi = double(I(allIndsFI));
        allIndsFJ = feat2(2,:,:) + (feat2(1,:,:)-1).*D.nRow;
        fj = double(I(allIndsFJ));

        % F*K x numTest logical vectors that assigns into 2^F
        signature = bsxfun(@gt, (fi - fj), [ferns.tau]); 
        signature = reshape(signature, [param.F param.K]);
        % divide into 2^F bins
        binAsg = sum(bsxfun(@times, signature, 2.^(0:param.F-1)'))...
                 + offset + 1;    
        outputs = [Rs(t).ferns.output];
        St = St + sum(outputs(:, binAsg), 2);
        
        % if VERBOSE & (t==1 | t==10)
        %     StSubNew = reshape(St, [2 D.nParts]);    
        %     subplot(1,2,(t==10)+1);
        %     visualize(I, reshape(initShapes(:, trial), [2 D.nParts]), StSubNew, D, ...
        %               feat1, feat2, Sgt);   
        %     title(sprintf('trial %d t=%d', trial, t));
        % end
    end
    % if VERBOSE
    %     err = sqrt(sum((Sgt - StSubNew).^2))./interocular;
    %     suptitle(sprintf('error %g', mean(err)));
    % end
    StAll(:, trial) = St;
end

ME = bsxfun(@rdivide, sqrt(sum(bsxfun(@minus,Sgt, reshape(initShapes,[2,D.nParts, ...
                   param.numTrials])).^2)), interocular);

StMedian = median(StAll, 2);
StMedianSub = reshape(median(StAll, 2), [2 D.nParts]);
fME = bsxfun(@rdivide, sqrt(sum((Sgt - StMedianSub).^2)), interocular);
fME = mean(fME);
if VERBOSE
fprintf('avg starting error=%.3g\t', mean(ME(:)));         
fprintf('final error=%.5g\tdelta=%.5g\n', fME, mean(ME(:))-fME);
end

