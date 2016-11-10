%TO DO: Make GUI
%TO DO: Remember that check_smear_cells only outputs an index corresponding
%to rows of cellPairs/cellImage that we want to keep
%TO DO: Maybe don't normalize mapped?
%TO DO Test a variety of TSNE conditions for plotting cells


caca;
folder_name = uigetdir('C:\Users\Jon\Desktop\');
file_list_rgb = dir(fullfile(folder_name,'*RGB*.tif'));
file_list_binary = dir(fullfile(folder_name,'*binary*.tif'));


for k = 1:length(file_list_rgb)
    file_names_rgb{k} = file_list_rgb(k).name;
    file_names_binary{k} = file_list_binary(k).name;
end

output = cell(length(file_names_rgb),1);

hWait = waitbar(0,'Segmenting blood...');
for k = 1:length(file_names_rgb)
   output{k} = segment_blood(fullfile(folder_name,file_names_rgb{k}),...
       fullfile(folder_name,file_names_binary{k}));
   
   waitbar(k/length(file_names_rgb));
end
close(hWait);

%%
for ii = 1:length(output)
    if ii == 1
        combined = output{ii}.cellPairs(:,3:12);
        combined = [myNorm(combined), ii*ones(length(output{ii}.cellPairs),1)];
    else
        combined = [combined; [myNorm(output{ii}.cellPairs(:,3:12)), ii*ones(length(output{ii}.cellPairs),1)]];
    end
    
end

groupIndex = combined(:,11);
combined(:,11) = [];

%%
sampleSize = 500;
for ii = 1:length(output)
    if ii == 1
        combinedRand = output{ii}.cellPairs(:,3:12);
        combinedRand = [myNorm(combinedRand), ii*ones(length(output{ii}.cellPairs),1)];
        combinedRand = datasample(combinedRand,sampleSize,'Replace',false);
    else
        tempCombined = [myNorm(output{ii}.cellPairs(:,3:12)), ii*ones(length(output{ii}.cellPairs),1)];
        combinedRand = [combinedRand; datasample(tempCombined,sampleSize,'Replace',false)];
    end
    
end

groupIndexRand = combinedRand(:,11);
combinedRand(:,11) = [];

%% Run TSNE
tsneSample = questdlg('All data or sample?','All data or sample','All data','Sample','Sample');

switch tsneSample
    case 'All data'
        data = combined;
    case 'Sample'
        data = combinedRand;
end

disp('Starting TSNE...');
include = 1:10;
% pcaIndex = pca(combined(:,include));

mappedX = fast_tsne(data(:,include),[],30,0.5);
mappedX(:,1) = myNorm(mappedX(:,1));
mappedX(:,2) = myNorm(mappedX(:,2));


%% Plot each
ss = 5000;
cellScaleFactor = 1;

switch tsneSample
    case 'All data'
        scatter_group = groupIndex;
    case 'Sample'
        scatter_group = groupIndexRand;
end

figure, gscatter(mappedX(:,1),mappedX(:,2),scatter_group,[],[],20);

myBloodSmearTSNEPlots([ones(length(data),1),ones(length(data),1),...
    data(:,include)],output{1}.cellPairLabels([1,2,include+2]),mappedX,ss);

for k = 1:length(output)
    
    myBloodSmearTSNEPlotsWithPics(output{k}.cellPairs(:,[1,2,include+2]),output{k}.cellImage,...
        mappedX(scatter_group == k,:),ss,cellScaleFactor,file_names_rgb{k});
    
%     myBloodSmearTSNEPlots(output{k}.cellPairs,output{k}.cellPairLabels,...
%         mappedX(groupIndex == k,:),ss);
end


