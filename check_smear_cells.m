function [idx] = check_smear_cells(cellPairs,cellImage,cellPairLabels,threshold)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% TO DO: Maybe have user input to choose specific files from a folder for
% the BloodSegmentation script? And if the output already exists, then you
% just append to it? That way if there is an error, you don't have to
% restart the whole thing. Should also then make sure that if there is an
% error, the sample is skipped, with an output saying what samples were
% skipped.
% TO DO: Add a Finish button to skip the rest of the cells.


%% Plot histograms and define outliers
nucCytMatch = cellPairs(:,1:2);
cellPairs(:,1:2) = [];

% if ~exist('threshold','var')
%     threshold = cellThresholds(cellPairs,cellPairLabels);
% end
% 
%     function threshold = cellThresholds(cellPairs,cellPairLabels)
%         for jj = 1:size(cellPairs,2)
%             figHandle = figure;
%             histogram(cellPairs(:,jj));
%             title(cellPairLabels{jj});
%             xl = xlim;
%             
%             [threshold{jj},~] = ginput(2);
%             threshold{jj} = sort(threshold{jj});
%             
%             
%             if isempty(threshold{jj}) == 1
%                 threshold{jj} = [0 xl(2)];
%             end
%             
%             if size(threshold{jj},1) == 1 & threshold{jj} < median(cellPairs(:,jj))
%                 threshold{jj}(2) = xl(2);
%                 
%             elseif size(threshold{jj},1) == 1 & threshold{jj} > median(cellPairs(:,jj))
%                 threshold{jj}(2) = threshold{jj};
%                 threshold{jj}(1) = 0;
%                 
%             end
%             
%             close(figHandle);
%         end
%     end
% 
% 
% %Get indices of cellPairs that fall outside thresholds
% for nn = 1:size(cellPairs,2)
%     outsideThresh(:,nn) = cellPairs(:,nn)<threshold{nn}(1) |...
%         cellPairs(:,nn)>threshold{nn}(2);
% end
% 
% threshTest = sum(outsideThresh,2);
% threshIndex = find(threshTest == 0);
% 
% cellImage = cellImage(cellPairs(threshIndex,1));
% cellPairs = cellPairs(threshIndex,:);
% 
% anotherRound = questdlg('Draw another round of thresholds?','Yes Or No','Yes','No','Yes');
% switch anotherRound
%     case 'Yes'
%         threshold = cellThresholds(cellPairs,cellPairLabels);
%     case 'No'
% end
% clear outsideThresh threshTest threshIndex
% for nn = 1:size(cellPairs,2)
%     outsideThresh(:,nn) = cellPairs(:,nn)<threshold{nn}(1) |...
%         cellPairs(:,nn)>threshold{nn}(2);
% end
% 
% newThresholds = threshold;
% threshTest = sum(outsideThresh,2);
% threshIndex = find(threshTest == 0);
% 
% cellImage = cellImage(cellPairs(threshIndex,1));
% cellPairs = cellPairs(threshIndex,:);
% 
% p = numSubplots(size(cellPairs,2));
% figure,
% for ii = 1:size(cellPairs,2)
%     subplot(p(1),p(2),ii)
%     histogram(cellPairs(:,ii));
%     yl = ylim;
%     hold on
%     plot([threshold{ii}(1),threshold{ii}(1)],[0,yl(2)],'-r');
%     plot([threshold{ii}(2),threshold{ii}(2)],[0,yl(2)],'-r');
%     title(cellPairLabels{ii});
% end
% 
% 




%% Sort cell images by size of image
totalPixels = cellfun(@(x) size(x,1)*size(x,2),cellImage);
[~, imageSortOrder] = sort(totalPixels);

% cellImageCopy = cellImage;
% cellPairsCopy = cellPairs;

cellImage = cellImage(imageSortOrder);
% cellPairs = cellPairs(imageSortOrder,:);


%%
screenSize = get(groot,'Screensize');
screen_border = 100;
squareSize = 110;
imageSize = 100;

%% Find number of squares that will fit into screen border, and round down.
xx = ((screenSize(3)-2*screen_border)*(screenSize(4)-2*screen_border))/(squareSize*squareSize);
xx = floor(xx);
p = numSubplots(xx); %Find optimal arrangment of squares
figure('OuterPosition',screenSize);

pageTotal = round(numel(cellImage)./xx);
pageText = uicontrol('Style', 'text',...
    'String', [num2str(1), ' of ',num2str(pageTotal)],...
    'Position', [screenSize(3)/2 screenSize(4)-140 130 30],...
    'FontSize',20);


leftPositions = screen_border/2:squareSize:(squareSize*(p(2)-1)+screen_border/2);
bottomPositions = screen_border/2:squareSize:(squareSize*(p(1)-1)+screen_border/2);
positions = combvec(leftPositions, bottomPositions);
positions = positions';

imageCountPerPage = [1,xx];

cell_toggle = cell(xx,1);
next_button = uicontrol('Style','pushbutton','String','Next',...
    'Position',[screenSize(3)/2-250 screenSize(4)-150 130 50],...
    'FontSize',20,...
    'Callback',{@pushbutton_callback,pageText});
cellIndex = zeros(numel(cellImage),1);

for k = 1:numel(cellImage)
    
    if k < imageCountPerPage(2)
        cell_toggle{mod(k,xx)} = uicontrol('Style','togglebutton','Callback',{@toggle_callback,cell_toggle});
        cell_toggle{mod(k,xx)}.Max = 1;
        cell_toggle{mod(k,xx)}.Min = 0;
        cell_toggle{mod(k,xx)}.CData = imresize(cellImage{k},[imageSize imageSize]);
        cell_toggle{mod(k,xx)}.Position = [positions(mod(k,xx),1),positions(mod(k,xx),2),imageSize,imageSize];
        
    elseif k == imageCountPerPage(2)
        cell_toggle{xx} = uicontrol('Style','togglebutton','Callback',{@toggle_callback,cell_toggle});
        cell_toggle{xx}.Max = 1;
        cell_toggle{xx}.Min = 0;
        cell_toggle{xx}.CData = imresize(cellImage{k},[imageSize imageSize]);
        cell_toggle{xx}.Position = [positions(xx,1),positions(xx,2),imageSize,imageSize];
        
        uiwait;
        for ii = 1:numel(cell_toggle)
            if ii == 1
                cellIndex(imageCountPerPage(1)) = cell_toggle{ii}.Value;
            else
                cellIndex(imageCountPerPage(1)+ii-1) = cell_toggle{ii}.Value;
            end
            
        end
        
        imageCountPerPage(1) = imageCountPerPage(2) + 1;
        imageCountPerPage(2) = imageCountPerPage(2) + xx;
    end
    

    
    
end

cellIndex = 1 - cellIndex; %This makes cells with a value of 1 the cells we want to keep
[~, sortMappingIndex] = sort(imageSortOrder); %Map the sorted index back to the original
idx = cellIndex(sortMappingIndex); %Put cellIndex in order of the original data


% cellImageNew = cellImageCopy(cellIndex == 1);
% cellPairsNew = cellPairsCopy(cellIndex == 1,:);

close(gcf);

    function pushbutton_callback(varargin)
        h = varargin{3};
        hString = h.String;
        tf = find(isspace(hString) == 1);
        currentPageNumSize = tf(1)-1;
        
        v = str2double(hString(1:currentPageNumSize));
        
        if numel(num2str(v+1)) > numel(num2str(v))
            hString = strcat('x',hString);
            hString(1:currentPageNumSize+1) = num2str(v+1);
        else
            hString(1:currentPageNumSize) = num2str(v+1);
        end
        
        h.String = hString;
        uiresume(gcbf);
    end

    function toggle_callback(varargin)
        h = varargin{1};
        value = h.Value;
        colorData = h.CData;
        if value == 1
            h.CData = immultiply(colorData,0.2);
        else
            h.CData = immultiply(colorData,5);
        end
        
    end

end

