function [meanShape, meanBox, alignedShapes] = computeMeanShape_new(shapes, D, Is)
%%%%%%%%%%%%%%%%%%%%
% Compute meanshape by aligning the line between the pupils to be
% horizontal wrt to the midpoint
% Input:
%  - shapes: 2 by nParts by N
%  - D: struct holding index of fiducial points and number of parts
%    D.leye, D.reye, D.nParts
%
%%%%%%%%%%%%%%%%%%%%
N = size(shapes, 3);

%for debugging only
DISPLAY = 0;

center = [(D.nCol-1)/2, (D.nRow-1)/2]; % center in (x,y), the midpt
% dist = (norm(shapes(:, D.leye, 1)'-center)); %dist from midpt
%                                                  %to each eyes

shapes = bsxfun(@minus, shapes, center');
alignedShapes = shapes;
% target = [center; center-[dist,0]];
for i = 1:N
    I = Is(:,:,i);
    pt = shapes(:,:,i);
    [angle, ratio] = getRotationAndScale(pt, D);
    % I = imrotate(I, angle*180/pi, 'bilinear');
    % I = imresize(I, ratio, 'bilinear');
    % rotate and scale pts
    pt = ([cos(angle) sin(angle); -sin(angle) cos(angle)]*pt).*ratio;
    [newR, newC, ~] = size(I);
    alignedShapes(:,:,i) = bsxfun(@plus, pt, [newC; newR]*0.5);

    % sfigure(3); imshow(I, []); hold on; axis image;
    % plot(alignedShapes(1,[1,2],i), alignedShapes(2,[1,2],i), 'y-');
    % for k = 1:length(D.connectedParts)
    %     plot(alignedShapes(1, D.connectedParts{k},i), ...
    %          alignedShapes(2, D.connectedParts{k},i),'g.-','MarkerSize',14, ...
    %          'LineWidth', 1);  
    % end

end

meanShape = mean(alignedShapes, 3);
midpt = (meanShape(:, 1) + meanShape(:, 2))*0.5;
center = midpt + [0; 11];
xmin = center(1)-(D.nCol-1)/2;
ymin = center(2)-(D.nRow-1)/2;
xmax = center(1)+(D.nCol-1)/2;
ymax = center(2)+(D.nRow-1)/2;

% xmin = min(meanShape(1,:)); xmax = max(meanShape(1,:));
% ymin = min(meanShape(2,:)); ymax = max(meanShape(2,:));
% padTop = 30; padWidth = 30;
% meanBox = [xmin, ymin, width, height];
% meanBox = round([xmin-padWidth, ymin-padTop, (xmax-xmin)+2*padWidth, (ymax-ymin)+(padWidth+padTop)]);
imgBox = [xmin, ymin, xmax-xmin, ymax-ymin];
% extract features from at least 20px within the boundary
meanBox = [xmin+20, ymin+20, xmax-xmin-40, ymax-ymin-40];

if DISPLAY
    hfig = sfigure(3); clf; axis ij equal; hold on;
    plot(midpt(1), midpt(2), 'm*');
    plot(cntr(1), cntr(2), 'g*');
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
    % saveTightFigure(hfig, 'results/images/meanshape.png');
end


%% mean shape before
% badMean = mean(shapes, 3);
% xminb = min(badMean(1,:)); xmaxb = max(badMean(1,:));
% yminb = min(badMean(2,:)); ymaxb = max(badMean(2,:));
% hfig = sfigure; axis([xminb-5 xmaxb+5 yminb-5 ymaxb+5]); axis ij equal; hold on;
% for i = 1:length(D.connectedParts)
%     plot(badMean(1, D.connectedParts{i}), ...
%          badMean(2, D.connectedParts{i}),'b.-','MarkerSize',14, ...
%          'LineWidth', 1);  
% end
% title('mean shape without alignment');
% % print(hfig, '-dpng', 'results/meanshapeBefore.png');
function [angle, ratio] = getRotationAndScale(pts, D)

leye = pts(:,D.leye);
reye = pts(:,D.reye);

dy = leye(2) - reye(2);
dx = leye(1) - reye(1);

angle = atan(dy/dx);

ratio = D.eyeDist./norm(leye-reye);%dist

if ratio > 3 || ratio < 0.5
    keyboard
end


