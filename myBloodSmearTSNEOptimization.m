%ydata = tsne(X, labels, no_dims, initial_dims, perplexity);
%mappedX = fast_tsne(X, initial_dims, perplexity, theta)

perRange = [10,50];
perStepSize = 10;
thetaRange = [0.1,0.9];
thetaStepSize = 0.2;
perSteps = (perRange(2)-perRange(1))/perStepSize + 1;
thetaSteps = (thetaRange(2) - thetaRange(1))/thetaStepSize + 1;

% p = numSubplots(perSteps.*thetaSteps);
mappedX = cell(perSteps,thetaSteps);

figure,
plotCount = 1;
thCount = 1;
pCount = 1;
for ii = perRange(1):perStepSize:perRange(2)
    for kk = thetaRange(1):thetaStepSize:thetaRange(2)
        mappedX{pCount,thCount} = fast_tsne(cellPairsNew,[],ii,kk);
        
        subplot(perSteps,thetaSteps,plotCount),
        scatterHandle{pCount,thCount} = scatter(mappedX{pCount,thCount}(:,1),mappedX{pCount,thCount}(:,2),20);
        title(['P = ',num2str(ii),' ','Th = ', num2str(kk)]);
        hold on
        thCount = thCount + 1;
        plotCount = plotCount + 1;
    end
    pCount = pCount + 1;
%     plotCount = plotCount + 1;
end

hold off;