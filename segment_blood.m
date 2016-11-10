function [output] = segment_blood(rgb_file,bw_file)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ~exist('rgb_file','var')
    [spinFileName, spinPathName] = uigetfile('/*RGB*.tif');
    [spinBWFileName, spinBWPathName] = uigetfile([spinPathName,'/*binary*.tif']);
    spin = imread(fullfile(spinPathName,spinFileName));
    spinBW = imread(fullfile(spinBWPathName, spinBWFileName));
    if size(spin,1) > 15000
        spin = spin(1:15000,:,:);
        spinBW = spinBW(1:15000,:,:);
    end
    
    if size(spin,2) > 15000
        spin = spin(:,1:15000,:);
        spinBW = spinBW(:,1:15000,:);
    end
    
else
    spin = imread(rgb_file);
    spinBW = imread(bw_file);
    
    if size(spin,1) > 15000
        spin = spin(1:15000,:,:);
        spinBW = spinBW(1:15000,:,:);
    end
    
    if size(spin,2) > 15000
        spin = spin(:,1:15000,:);
        spinBW = spinBW(:,1:15000,:);
    end
    
    
    
end

%% Read in RGB and binary image files (from ImageJ)

% spin = imresize(spin,scaleFactor);

spinBW = imbinarize(spinBW);
spinBW = imcomplement(spinBW);
% spinBW = imresize(spinBW,scaleFactor);

spin = bsxfun(@times, spin, cast(spinBW,'like',spin));

clear spinBW
spinSize = size(spin);
spinSize(3) = [];

% Convert to HSV
spinHSV = rgb2hsv(spin);
spinH = spinHSV(:,:,1);
spinS = spinHSV(:,:,2);
% spinRed = spin(:,:,1);
% spinGreen = spin(:,:,2);
% spinBlue = spin(:,:,3);

%% K-means cluster the colors in RGB space

[ab, blackPixels, cluster_center, clusterSort,clusteridx] = kmeans_blood(spin,3);

    function [ab, blackPixels,cluster_center, clusterSort,clusteridx] = kmeans_blood(spin_image,nColors)
        ab = double(spin_image);
        clear spin_image;
        nrows = size(ab,1);
        ncols = size(ab,2);
        ab = reshape(ab,nrows*ncols,3);
        
        blackPixels = ab(:,1) == 0 & ab(:,2) == 0 & ab(:,3) == 0;
        nonBlackPixels = ~blackPixels;
        colorPixels = ab(nonBlackPixels == 1,:);
        pixelSample = datasample(colorPixels,round(0.01*length(colorPixels)),'Replace',false);
        [clusteridx, cluster_center] = kmeans(pixelSample, nColors,'distance','sqEuclidean',...
            'Replicates',3,'Start','cluster');
        
        [~, clusterSort] = sortrows(cluster_center,1);
    end



%% Nucleus is cluster center with min Red (clusterSort(1))
nucCluster = clusterSort(1);

nucStDevR = std(ab(clusteridx == nucCluster,1));
nucStDevG = std(ab(clusteridx == nucCluster,2));
nucStDevB = std(ab(clusteridx == nucCluster,3));

distThreshold = mean([nucStDevR, nucStDevG, nucStDevB]);

nucMeanR = cluster_center(nucCluster,1);
nucMeanG = cluster_center(nucCluster,2);
nucMeanB = cluster_center(nucCluster,3);

%% Find parameters that fit nucleus

nucImageBW = nucleus_param(spin,ab,blackPixels,[nucMeanR nucMeanG nucMeanB],spinSize,distThreshold);
clear ab;
clear blackPixels;

nucBoundaries = edge(nucImageBW);
% nucBoundaries = bwboundaries(nucImageBW,8,'noholes');
% 
% nucBoundariesMask = zeros(spinSize);
% for kk = 1:length(nucBoundaries)
%     boundary = nucBoundaries{kk};
%     nucBoundariesMask(boundary(:,2),boundary(:,1)) = 1;
% end


    function [nucImageBW] = nucleus_param(spin_image,ab, blackPixels,nucMeans, spinSize, distThreshold)
        nucPixels = zeros(size(ab,1),1);
        for k = 1:length(ab)
            if blackPixels(k) == 1
                nucPixels(k) = 0;
            else
                colorDist = norm(nucMeans-[ab(k,1),ab(k,2),ab(k,3)]);
                if colorDist < distThreshold
                    nucPixels(k) = 1;
                end
            end
            
        end
        clear blackPixels;
        clear ab;
        
        nucPixels = logical(nucPixels);
        % nucPixels = ~nucPixels;
        nucPixels2D = reshape(nucPixels,spinSize(1),spinSize(2));
        clear nucPixels;
        nucMask = cat(3,nucPixels2D,nucPixels2D,nucPixels2D);
        clear nucPixels2D;
        
        nucImage = bsxfun(@times, spin_image, cast(nucMask,'like',spin_image));
        clear spin_image;
        clear nucMask;
        % figure, imshow(nucImage);
        
        nucImageGray = rgb2gray(nucImage);
        clear nucImage;
        bwEdges = edge(nucImageGray,'sobel');
        
        % figure, imshow(bwEdges);
        
        %Dilate and fill nucleus, then erode the dilation
        bwEdgesDil = imdilate(bwEdges,strel('disk',1));
        clear bwEdges;
        bwEdgesFill = imfill(bwEdgesDil,'holes');
        nucImageBW = imerode(bwEdgesFill,strel('disk',3));
        % figure, imshow(nucImageBW),title('nucImageBW');
        
    end
%% Get properties of nucleus and cytoplasm

nucObjects = bwconncomp(nucImageBW);
nucProps = regionprops(nucObjects,'Area','Centroid','PixelIdxList','BoundingBox','Eccentricity');
nucPropsH = regionprops(nucObjects,spinH,'MeanIntensity');
nucPropsS = regionprops(nucObjects,spinS,'MeanIntensity');
% nucPropsR = regionprops(nucObjects,spinRed,'MeanIntensity');
% nucPropsG = regionprops(nucObjects,spinGreen,'MeanIntensity');
% nucPropsB = regionprops(nucObjects,spinBlue,'MeanIntensity');


%Cytoplasm connected components
cytImage = bsxfun(@times, spin, cast(imcomplement(nucImageBW), 'like', spin));
clear nucImageBW;

cytImageHSV = rgb2hsv(cytImage);
cytImageHSV = uint8(cytImageHSV.*256);

cytImageNaNV = cytImageHSV(:,:,3);
cytImageNaNV(cytImageNaNV == 0) = NaN;
cytImageNaNH = cytImageHSV(:,:,1);
cytImageNaNH(cytImageNaNV == 0) = NaN;
cytImageNaNS = cytImageHSV(:,:,2);
cytImageNaNS(cytImageNaNV == 0) = NaN;

clear cytImageHSV;

% cytImageNaNR = cytImage(:,:,1);
% cytImageNaNR(cytImageNaNR == 0) = NaN;
% cytImageNaNG = cytImage(:,:,2);
% cytImageNaNG(cytImageNaNG == 0) = NaN;
% cytImageNaNB = cytImage(:,:,3);
% cytImageNaNB(cytImageNaNB == 0) = NaN;

% figure, imshow(cytImage);

cytImage = rgb2gray(cytImage);

cytImageBW = imbinarize(cytImage,0.1);
clear cytImage;
cytImageBW = imfill(cytImageBW,'holes');


% cytBoundaries = bwboundaries(cytImageBW,8,'noholes');
% cytBoundariesMask = zeros(spinSize);
% for kk = 1:length(cytBoundaries)
%     boundary = cytBoundaries{kk};
%     cytBoundariesMask(boundary(:,2),boundary(:,1)) = 1;
% end
cytBoundaries = edge(cytImageBW);

cytObjects = bwconncomp(cytImageBW);

%Cytoplasm properties
cytProps = regionprops(cytObjects,'Area','Centroid','PixelIdxList','BoundingBox','Eccentricity');

cytPropsH = myregionprops(cytObjects,cytImageNaNH,'MeanIntensity');
cytPropsS = myregionprops(cytObjects,cytImageNaNS,'MeanIntensity');

% cytPropsR = myregionprops(cytObjects,cytImageNaNR,'MeanIntensity');
% cytPropsG = myregionprops(cytObjects,cytImageNaNG,'MeanIntensity');
% cytPropsB = myregionprops(cytObjects,cytImageNaNB,'MeanIntensity');

%% Get cells with only one nuclei

[cellPairs,oneNuc, oneNucBoundBox] = one_nuc(nucProps,cytProps);

    function [cellPairs,oneNuc, oneNucBoundBox] = one_nuc(nucProps,cytProps)
        boundBox = {cytProps.BoundingBox}';
        clear cytProps;
        boundBox = cell2mat(boundBox);
        
        % boundBox(areaIndex,:) = [];
        boundBox = double(boundBox);
        
        nucCentroid = {nucProps.Centroid}';
        clear nucProps;
        nucCentroid = cell2mat(nucCentroid);
        
        pointsToTest = nucCentroid;
        
        x1 = round(boundBox(:,1));
        x2 = round(boundBox(:,1)+boundBox(:,3));
        y1 = round(boundBox(:,2));
        y2 = round(boundBox(:,2)+boundBox(:,4));
        
        oneNuc = zeros(length(boundBox),1);
        for jj = 1:length(boundBox)
            
            node = [x1(jj) y1(jj); x2(jj) y1(jj); x2(jj) y2(jj); x1(jj) y2(jj); x1(jj) y1(jj)];
            test = inpoly(pointsToTest,node);
            
            if sum(test) == 1
                oneNuc(jj) = 1;
                nucIndex(jj) = find(test == 1);
            end
            
        end
        
        nucIndex(nucIndex == 0) = [];
        nucIndex = nucIndex';
        count = 1:length(nucIndex);
        count = count';
        
        cellPairs = [count,nucIndex];
        oneNucBoundBox = boundBox(oneNuc == 1,:);
    end

%% Cell images
cellImage = get_cell_images(oneNucBoundBox,spin,spinSize,nucBoundaries,cytBoundaries);

    function cellImage = get_cell_images(oneNucBoundBox,spin,spinSize,nucBoundaries, cytBoundaries)
%         subs = numSubplots(length(oneNucBoundBox));
        pixelPad = 5;
        
        cellImage = cell(length(oneNucBoundBox),1);
        
        % figure,
        for r = 1:length(oneNucBoundBox)
            %     subplot(subs(1),subs(2),ii),
            x = round(oneNucBoundBox(r,1)-pixelPad:(oneNucBoundBox(r,1)+oneNucBoundBox(r,3)+pixelPad));
            x(x<=0) = 1;
            x(x>spinSize(2))=spinSize(2);
            y = round(oneNucBoundBox(r,2)-pixelPad:(oneNucBoundBox(r,2)+oneNucBoundBox(r,4)+pixelPad));
            y(y<=0) = 1;
            y(y>spinSize(1)) = spinSize(1);
            cellImageRed = spin(y,x,1);
            cellImageGreen = spin(y,x,2);
            cellImageBlue = spin(y,x,3);
            nucMask = nucBoundaries(y,x);
            cytMask = cytBoundaries(y,x);
            cellImageRed(cytMask) = 255;
            cellImageGreen(cytMask) = 0;
            cellImageBlue(cytMask) = 0;
            cellImageRed(nucMask) = 255;
            cellImageGreen(nucMask) = 255;
            cellImageBlue(nucMask) = 0;
            
            cellImage{r} = cat(3,cellImageRed,cellImageGreen,cellImageBlue);
            
            %     imshow(myImage(y,x,:),[]);
        end
    end

%% Properties to measure
[cellPairs, cellPairLabels] = cell_properties(oneNuc,cellPairs,nucProps,cytProps,nucPropsH,nucPropsS,...
    cytPropsH,cytPropsS);

    function [cellPairs, cellPairLabels] = cell_properties(oneNuc,cellPairs,nucProps,cytProps,...
            nucPropsH,nucPropsS,cytPropsH,cytPropsS)
        
        %Nucleus
        nucArea = [nucProps.Area]';
        nucCentroid = {nucProps.Centroid}';
        nucCentroid = cell2mat(nucCentroid);
        nucEccent = [nucProps.Eccentricity]';
        nucH = [nucPropsH.MeanIntensity]';
        nucS = [nucPropsS.MeanIntensity]';
        
%         nucR = [nucPropsR.MeanIntensity]';
%         nucG = [nucPropsG.MeanIntensity]';
%         nucB = [nucPropsB.MeanIntensity]';
        
        combinedNuc = [nucArea, nucEccent, nucH, nucS, nucCentroid];
        combinedNuc = combinedNuc(cellPairs(:,2),:);
        
        %Cytoplasm
        cytArea = [cytProps.Area]';
        cytCentroid = {cytProps.Centroid}';
        cytCentroid = cell2mat(cytCentroid);
        cytEccent = [cytProps.Eccentricity]';
        cytH = [cytPropsH.MeanIntensity]';
        cytS = [cytPropsS.MeanIntensity]';
        
%         cytR = [cytPropsR.MeanIntensity]';
%         cytG = [cytPropsG.MeanIntensity]';
%         cytB = [cytPropsB.MeanIntensity]';
        
        combinedCyt = [cytArea,cytEccent, cytH, cytS,cytCentroid];
        
        combinedCyt = combinedCyt(oneNuc == 1,:);
        
        cellPairLabels = {'Cyto index','NucIndex','Cyt Area','Nuc Area',...
            'Nuc/Cyt Area Ratio','Cyt Mean Color H', 'Cyt Mean Color S',...
            'Nuc Mean Color H', 'Nuc Mean Color S',...
            'Centroid distance', 'Cyt Eccentricity', 'Nuc Eccentricity'};
        %1 Cyto index
        
        %2 Nuc Index
        
        %3 cytArea
        cellPairs(:,3) = combinedCyt(:,1);
        %4 nucArea
        cellPairs(:,4) = combinedNuc(:,1);
        %5 nucArea/cytArea ratio
        cellPairs(:,5) = cellPairs(:,4)./cellPairs(:,3);
        %6 cytH mean color H
        cellPairs(:,6) = combinedCyt(:,3);
        %7 cytS mean color S
        cellPairs(:,7) = combinedCyt(:,4);
        %8 nucH mean color H
        cellPairs(:,8) = combinedNuc(:,3);
        %9 nucS mean color S
        cellPairs(:,9) = combinedNuc(:,4);
        %10 distance between centroids
        c1 = bsxfun(@minus,combinedCyt(:,5),combinedNuc(:,5));
        c2 = bsxfun(@minus,combinedCyt(:,6),combinedNuc(:,6));
        distances = sqrt(c1.^2+c2.^2);
        cellPairs(:,10) = distances;
        %11 cytEccent
        cellPairs(:,11) = combinedCyt(:,2);
        %12 nucEccent
        cellPairs(:,12) = combinedNuc(:,2);
    end

%% Exclusion properties
% cellPairs(cellPairs(:,5)>1,:) = []; %'Cells' that have nuc bigger than cyt disregarded
% cellPairs(cellPairs(:,5)<0.1,:) = []; % Cells that have very small nuclei disregarded
% cellPairs(cellPairs(:,3)>1e4,:) = [];

% p = numSubplots(size(cellPairs,2)-2);
% figure,
% for ii = 1:size(cellPairs,2)-2
%     subplot(p(1),p(2),ii)
%     histogram(cellPairs(:,ii+2));
%     title(cellPairLabels{ii+2});
% end

% %% Run TSNE
% pcaIndex = pca(cellPairs(:,3:11));
% 
% mappedX = fast_tsne(cellPairs(:,3:11),[],30,0.1);
% mappedX(:,1) = myNorm(mappedX(:,1));
% mappedX(:,2) = myNorm(mappedX(:,2));
% % aspectRatio = max(mappedX(:,1))/max(mappedX(:,2));
% 
% % Plot
% ss = 5000;
% cellScaleFactor = 0.5;
% myBloodSmearTSNEPlotsWithPics(cellPairs,cellImage,mappedX,ss,cellScaleFactor);
% 
% myBloodSmearTSNEPlots(cellPairs,cellPairLabels,mappedX,ss);

% %% Cluster cell populations
% clusterIndex = clusterdata(mappedX,15);
% % clusterIndex = clusterdata(cellPairs(:,3:11),3);
% colorsForClustering = distinguishable_colors(max(clusterIndex(:)));
% figure,
% scatter(mappedX(:,1),mappedX(:,2),15,colorsForClustering(clusterIndex,:),'filled');
% % scatter(pcaIndex(:,1),pcaIndex(:,2),15,colorsForClustering(clusterIndex,:),'filled');

%% Generate output structure
output.cellImage = cellImage;
output.cellPairs = cellPairs;
output.cellPairLabels = cellPairLabels;
% output.pcaIndex = pcaIndex;
% output.mappedX = mappedX;


end

