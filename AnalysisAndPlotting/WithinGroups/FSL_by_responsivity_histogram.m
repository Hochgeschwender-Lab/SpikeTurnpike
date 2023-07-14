function FSL_by_responsivity_histogram(all_data,trialTagsLabels,trialType)
%AUROC_HISTOGRAM Summary of this function goes here
%   Detailed explanation goes here

groupNames = fieldnames(all_data);

trialTypeInd = find(strcmp(trialTagsLabels,trialType));

FSL_vec = [];
cellTypesVec = {};
responsivityVec = {};

for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};
    mouseNames = fieldnames(all_data.(groupName));

    for mouseNum = 1:length(mouseNames)
        mouseName = mouseNames{mouseNum};
        cellIDs = fieldnames(all_data.(groupName).(mouseName));

        for cellID_num = 1:length(cellIDs)
            cellID = cellIDs{cellID_num};

            FSL_vec(end+1,1) = all_data.(groupName).(mouseName).(cellID).FirstSpikeLatency(trialTypeInd);

            cellTypesVec{end+1,1} = all_data.(groupName).(mouseName).(cellID).Cell_Type;

            responsivityNum = all_data.(groupName).(mouseName).(cellID).StimResponsivity;
            if responsivityNum == 1
                responsivityVec{end+1,1} = '+';
            elseif responsivityNum == 0
                responsivityVec{end+1,1} = 'NR';
            else
                responsivityVec{end+1,1} = '-';
            end
        end
    end
end

%% Plotting
cell_types = unique(cellTypesVec);

figure;
T = tiledlayout('flow');
for cellTypeNum = 1:length(cell_types)
    FSL_vec_cellType = FSL_vec(strcmp(cellTypesVec, cell_types{cellTypeNum}));
    responsivityVec_cellType = responsivityVec(strcmp(cellTypesVec, cell_types{cellTypeNum}));

    t = nexttile(T);
    title(t, cell_types{cellTypeNum});
    hold on
    h1 = histogram(FSL_vec_cellType(strcmp(responsivityVec_cellType,'+')), 'BinWidth',5);
    h2 = histogram(FSL_vec_cellType(strcmp(responsivityVec_cellType,'NR')), 'BinWidth',5);
    % h3 = histogram(AUROC_vec_cellType(strcmp(responsivityVec_cellType,'-')), 20, 'BinLimits', [0 1], 'FaceColor','red', 'FaceAlpha',0.5);
   
    % [hc,he]=histcounts(AUROC_vec_cellType(strcmp(responsivityVec_cellType,'+')), 0:0.04:1);
    % hc(end+1) = 0;
    % stairs(he, hc, 'Color', 'green', 'LineWidth',4);
    % 
    % [hc,he]=histcounts(AUROC_vec_cellType(strcmp(responsivityVec_cellType,'NR')), 0:0.04:1);
    % hc(end+1) = 0;
    % stairs(he, hc, 'Color', 'yellow', 'LineWidth',4);
    % 
    % [hc,he]=histcounts(AUROC_vec_cellType(strcmp(responsivityVec_cellType,'-')), 0:0.04:1);
    % hc(end+1) = 0;
    % stairs(he, hc, 'Color', 'red', 'LineWidth',4);
    
    hold off

    xlabel('First Spike Latency (ms)');
    ylabel('Unit Count');
end

end

