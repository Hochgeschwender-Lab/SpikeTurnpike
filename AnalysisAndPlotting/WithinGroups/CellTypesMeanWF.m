function CellTypesMeanWF(all_data,cell_types,cell_types_to_plot)
%CELLTYPESMEANWF Summary of this function goes here
%   Detailed explanation goes here

groupNames = fieldnames(all_data);
colors = brewermap(length(cell_types),'Set2');

CellTypeInds_vec = []; % vector of all cell types
waveforms = []; % [samples x units]
groupsLabels = {};
mouseLabels = {};
cellID_Labels = {};

for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};
    mouseNames = fieldnames(all_data.(groupName));

    for mouseNum = 1:length(mouseNames)
        mouseName = mouseNames{mouseNum};
        cellIDs = fieldnames(all_data.(groupName).(mouseName));

        for cellID_num = 1:length(cellIDs)
            cellID = cellIDs{cellID_num};
            cellType = all_data.(groupName).(mouseName).(cellID).Cell_Type;
            if strcmp(cellType,'TAN') % delete this!
                cellType = 'MSN';
            end
            isSingleUnit = all_data.(groupName).(mouseName).(cellID).IsSingleUnit;
            if any(strcmp(cell_types_to_plot,cellType)) && isSingleUnit
                cellTypeInd = find(strcmp(cell_types,cellType));
                CellTypeInds_vec = [CellTypeInds_vec; cellTypeInd];

                waveform = all_data.(groupName).(mouseName).(cellID).Normalized_Template_Waveform;
                waveforms = [waveforms, waveform];

                groupsLabels{end+1,1} = groupName;
                mouseLabels{end+1,1} = mouseName;
                cellID_Labels{end+1,1} = cellID;
            end
        end
    end
end

%% Plotting

%fig = figure();
% for unitInd = 1:size(waveforms,2)
%      plot(waveforms(:,unitInd), 'Color',colors(CellTypeInds_vec(unitInd),:));
% end
figure();
T = tiledlayout(length(cell_types_to_plot)+1,1);
meanWFs = [];
colorInds = [];
for cellTypeInd = 1:length(cell_types)
    if any(strcmp(cell_types_to_plot, cell_types{cellTypeInd}))
        nexttile(T);
    
        inds_thisCellType = find(CellTypeInds_vec==cellTypeInd);
        colorInds(end+1) = cellTypeInd;
    
        waveforms_thisCellType = waveforms(:,inds_thisCellType);
        groupsLabels_thisCellType = groupsLabels(inds_thisCellType);
        mouseLabels_thisCellType = mouseLabels(inds_thisCellType);
        cellID_Labels_thisCellType = cellID_Labels(inds_thisCellType);
    
        meanWF = mean(waveforms_thisCellType,2);
        meanWFs = [meanWFs meanWF];
        hold on
        for unitInd = 1:size(waveforms_thisCellType,2)
            p = plot(waveforms_thisCellType(:,unitInd), 'Color',colors(cellTypeInd,:));
    
            % Data tips so sus units can be identified
            p.DataTipTemplate.DataTipRows(1:2) = [dataTipTextRow("Rec:", repmat(mouseLabels_thisCellType(unitInd),size(waveforms_thisCellType,1))),...
                dataTipTextRow("Cell:", repmat(cellID_Labels_thisCellType(unitInd),size(waveforms_thisCellType,1)))];
            p.DataTipTemplate.set('Interpreter','none');
        end
        plot(meanWF, 'Color','k');
        hold off
    end
end

nexttile(T);
hold on
for meanWF_ind = 1:size(meanWFs,2)
    plot(meanWFs(:,meanWF_ind), 'Color',colors(colorInds(meanWF_ind),:));
end
hold off

legend(cell_types_to_plot);

end

