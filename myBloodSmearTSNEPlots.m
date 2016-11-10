function [scatterHandle] = myBloodSmearTSNEPlots(cellPairs,cellPairLabels,mappedX,ss,title_name)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% TSNE plots for each metric
p = numSubplots(size(cellPairs,2)-2);
figure,
for ii = 1:size(cellPairs,2)-2
    subplot(p(1),p(2),ii)
    scatterHandle(ii) = scatter(ss.*mappedX(:,1),ss.*mappedX(:,2),10,cellPairs(:,ii+2),'filled');
%     axis('equal')
    title(cellPairLabels{ii+2});
    set(gca,'XTick',[],'YTick',[]);
    hold on
end
hold off
colormap('jet');
if exist('title_name','var')
    title(title_name,'Interpreter','none');
end

end

