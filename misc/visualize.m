function visualize(I, Sold, St, D,f1, f2, Sgt)
    if length(size(I))==3
        I = reshape(I, [D.nRow D.nCol 3]);    
    else    
        I = reshape(I, [D.nRow D.nCol]);
    end
    imshow(I, []); hold on;
    if nargin > 4
        if ~isempty(f1)
            colors = {'r','c','y','w','m'};
            nSets = size(f1,2)/5;
            for j = 1:nSets
                start = (j-1)*5;
                for i=1:5
                    plot(f1(1,start+i),f1(2,start+i), [colors{i},'+'], 'MarkerSize', 7);
                    plot(f2(1,start+i),f2(2,start+i), [colors{i},'o'], 'MarkerSize', 7);
                end
            end
        end
    end
    for i=1:length(D.connectedParts)
        plot(Sold(1, D.connectedParts{i}), ...
                 Sold(2, D.connectedParts{i}),'r.-','MarkerSize',8);
        plot(St(1, D.connectedParts{i}), ...
                 St(2, D.connectedParts{i}),'g.-','MarkerSize',8);
            % plot(Sold(1, D.connectedParts{i}), ...
            %      Sold(2, D.connectedParts{i}),'k-','MarkerSize',8, ...
            %      'LineWidth', 2.2);
    end
    if nargin > 6
        for i = 1:length(D.connectedParts)
        % plot(Sgt(1, D.connectedParts2{i}), ...
        %          Sgt(2, D.connectedParts2{i}),'c.','MarkerSize',14);  
        plot(Sgt(1, D.connectedParts{i}), ...
             Sgt(2, D.connectedParts{i}),'b.--','MarkerSize',7, ...
             'LineWidth', 0.5);  
        % plot(Sgt(1, D.connectedParts{i}), ...
        %      Sgt(2, D.connectedParts{i}),'b.-','MarkerSize',8, ...
        %      'LineWidth', 1);  

        end
    end
    % leyebrow = [1,2,4,6,8,7,5,3];
    % leyeTop = [9 10 11];
    % leye = [14 16 12 13 15];
    % reyebrow = [17 18 20 22 24 23 21 19];
    % reyeTop = [27 26 25];
    % reye = [28 29 30 32 31];
    % nose1 = [33 34 35]; 
    % nose2 = [34 36 38 40];
    % nose3 = [37 38 39];
    % mouthout = [41 42 44 46 48 54 52 50];
    % mouthin = [43 45 47 53 51 49];
