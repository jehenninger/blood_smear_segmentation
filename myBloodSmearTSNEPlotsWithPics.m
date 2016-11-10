function myBloodSmearTSNEPlotsWithPics(cellPairs,cellImage,mappedX,ss,cellScaleFactor,title_name)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

figure,
% scatter(ss*mappedX(:,1),ss*mappedX(:,2));
fudgeFactor = 0;
hold on
for ii = 1:size(cellPairs,1)
%     cellPairImage = cellImage{cellPairs(ii,1)};
    cellPairImage = cellImage{ii};
    scaledCellImage = imresize(cellPairImage,[size(cellPairImage,1)*cellScaleFactor,...
        size(cellPairImage,2)*cellScaleFactor]);
    image(ss*mappedX(ii,1)-fudgeFactor,ss*mappedX(ii,2)-fudgeFactor,scaledCellImage);
end
hold off
% xlim([0-0.01*ss ss+0.01*ss]);
% ylim([0-0.01*ss ss+0.01*ss]);
axis('equal');

if exist('title_name','var')
    title(title_name,'Interpreter','none');
end

end

