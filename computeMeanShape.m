function [meanShape, meanBox, alignedShapes] = computeMeanShape(shapes, D, Is)
%%%%%%%%%%%%%%%%%%%%
% Compute meanshape by aligning all shapes by having the origin at
% the midpoint between the two pupils, scaling the distance between
% that midpoint and the mouth to have a fixed constant, and by
% transforming that that line is vertical.
% Input:
%  - shapes: 2 by nParts by N
%  - D: struct holding index of fiducial points and number of parts
%    D.leye, D.reye, D.mouth, D.nParts
%
%%%%%%%%%%%%%%%%%%%%
N = size(shapes, 3);

%for debugging only
DISPLAY = 0;
if DISPLAY, Is = reshape(Is, [250*250, N]);, end

midpt = (shapes(:, D.leye, :) + shapes(:, D.reye, :))./2;
origin = mean(midpt, 3)'; 
mid2mouth = sqrt(sum((midpt - shapes(:, D.mouth, :)).^2));
scale = round(mean(mid2mouth));

if DISPLAY, fprintf('scale s.t. the distance from midpoint to mouth is %g\n', scale);end
alignedShapes = shapes;
target = [origin; origin+[0 scale]];
for i = 1:N
    T = bestfit_nonreflective_similarity(target,...
                                         [midpt(:,1,i)'; shapes(:, D.mouth, i)']);
    alignedShapes(:,:,i) = ([alignedShapes(:,:,i)' ones(D.nParts, 1)]*T)';
end

meanShape = mean(alignedShapes, 3);
xmin = min(meanShape(1,:)); xmax = max(meanShape(1,:));
ymin = min(meanShape(2,:)); ymax = max(meanShape(2,:));
padTop = 30; padWidth = 20;
meanBox = round([xmin-padWidth, ymin-padTop, (xmax-xmin)+2*padWidth, (ymax-ymin)+(padWidth+padTop)]);

if DISPLAY
    hfig = sfigure; axis([xmin-padWidth-5 xmax+padWidth+5 ymin-padTop-5 ymax+padWidth+5]); axis ij image; hold on;
    for i = 1:length(D.connectedParts)
        plot(meanShape(1, D.connectedParts{i}), ...
             meanShape(2, D.connectedParts{i}),'b.-','MarkerSize',14, ...
             'LineWidth', 1);  
    end
    rectangle('Position', [xmin, ymin, xmax-xmin, ymax-ymin], 'EdgeColor', ...
              'g', 'LineStyle', '--');    
    rectangle('Position', meanBox, 'EdgeColor', 'g', 'LineWidth', 3);    
    title('mean shape after alignment');    
    pinds = randsample(meanBox(4)*meanBox(3), 400); 
    ri = rem(pinds-1, meanBox(4))+1 + meanBox(2); 
    ci = (pinds-ri)/meanBox(4) +1 + meanBox(1);
    plot(ci, ri, 'c.');
    % print(hfig, '-dpng', 'results/meanshape.png');
end

%%%%%%%Draw on %%%%%%%%%%%%%
centerMean = mean(meanShape, 2)';

bboxW = floor(D.nRow/2.2);
imgOffset = (D.nRow-bboxW)/2;%from left corner of the image to the
                             %left corner of the bbox
if DISPLAY
    imgInd = randi(N, 9,1);
    sfigure;
    for i = 1:3
        vis = imgInd(i);
        subplot(1,3,i);visualize(Is(:,vis), shapes(:,:,vis), meanShape , ...
                                 D, [],[]);
        rectangle('Position', [imgOffset, imgOffset, bboxW, bboxW], 'EdgeColor', 'r',...
                  'LineStyle', '--');    
        rectangle('Position', meanBox, 'EdgeColor', 'g');    
        % pind = randsample(meanBox(4)*meanBox(3), 400); 
        % ri = rem(pind-1, meanBox(4))+1+ meanBoxAdjusted(2);
        % ci = (pind-ri)/meanBox(4)+1+ meanBoxAdjusted(1);
        % plot(ci, ri, 'm.');
    end
end
%% mean shape before
% badMean = mean(shapes, 3);
% xminb = min(badMean(1,:)); xmaxb = max(badMean(1,:));
% yminb = min(badMean(2,:)); ymaxb = max(badMean(2,:));
% pad = 20;
% hfig = sfigure; axis([xminb-padWidth-5 xmaxb+padWidth+5 yminb-padTop-5 ymaxb+padWidth+5]); axis ij image; hold on;
% for i = 1:length(D.connectedParts)
%     plot(badMean(1, D.connectedParts{i}), ...
%          badMean(2, D.connectedParts{i}),'b.-','MarkerSize',14, ...
%          'LineWidth', 1);  
% end
% rectangle('Position', [xminb, yminb, xmaxb-xminb, ymaxb-yminb], 'EdgeColor', ...
%           'g', 'LineStyle', '--');    
% badMeanBox = [xminb-padWidth, yminb-padTop, (xmaxb-xminb)+2*padWidth, (ymaxb-yminb)+(padWidth+padTop)];
% rectangle('Position', badMeanBox, 'EdgeColor', 'g', 'LineWidth', 3);    
% title('mean shape without alignment');
% print(hfig, '-dpng', 'results/meanshapeBefore.png');
