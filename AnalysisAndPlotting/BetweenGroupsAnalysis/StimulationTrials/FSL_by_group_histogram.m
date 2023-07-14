function FSL_by_group_histogram(all_data,trialTagsLabels,trialType)
% Histogram of only positively-modulated units' first spike latencies,
% overlayed between groups.

groupNames = fieldnames(all_data);

trialTypeInd = find(strcmp(trialTagsLabels,trialType));

FSL_vec = [];
cellTypesVec = {};
groupVec = {};

for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};
    mouseNames = fieldnames(all_data.(groupName));

    for mouseNum = 1:length(mouseNames)
        mouseName = mouseNames{mouseNum};
        cellIDs = fieldnames(all_data.(groupName).(mouseName));

        for cellID_num = 1:length(cellIDs)
            cellID = cellIDs{cellID_num};

            isSingleUnit = all_data.(groupName).(mouseName).(cellID).IsSingleUnit;
            responsivityNum = all_data.(groupName).(mouseName).(cellID).StimResponsivity;
            %if (responsivityNum == 1) && isSingleUnit % unit is positively modulated by stim
            if isSingleUnit && (responsivityNum == 1)
                groupVec{end+1,1} = groupName;
                FSL_vec(end+1,1) = all_data.(groupName).(mouseName).(cellID).FirstSpikeLatency(trialTypeInd);
                cellTypesVec{end+1,1} = all_data.(groupName).(mouseName).(cellID).Cell_Type;
            end
        end
    end
end

%% Plotting
cell_types = unique(cellTypesVec);

% figure;
% T = tiledlayout('flow');
% for cellTypeNum = 1:length(cell_types)
%     FSL_vec_cellType = FSL_vec(strcmp(cellTypesVec, cell_types{cellTypeNum}));
%     groupVec_cellType = groupVec(strcmp(cellTypesVec, cell_types{cellTypeNum}));
% 
%     t = nexttile(T);
%     title(t, cell_types{cellTypeNum});
%     hold on
%     h1 = histogram(FSL_vec_cellType(strcmp(groupVec_cellType,groupNames{1})), 'BinWidth',3, 'FaceColor','red');
%     h2 = histogram(FSL_vec_cellType(strcmp(groupVec_cellType,groupNames{2})), 'BinWidth',3, 'FaceColor','blue');
%     % h3 = histogram(AUROC_vec_cellType(strcmp(responsivityVec_cellType,'-')), 20, 'BinLimits', [0 1], 'FaceColor','red', 'FaceAlpha',0.5);
% 
%     % [hc,he]=histcounts(AUROC_vec_cellType(strcmp(responsivityVec_cellType,'+')), 0:0.04:1);
%     % hc(end+1) = 0;
%     % stairs(he, hc, 'Color', 'green', 'LineWidth',4);
%     % 
%     % [hc,he]=histcounts(AUROC_vec_cellType(strcmp(responsivityVec_cellType,'NR')), 0:0.04:1);
%     % hc(end+1) = 0;
%     % stairs(he, hc, 'Color', 'yellow', 'LineWidth',4);
%     % 
%     % [hc,he]=histcounts(AUROC_vec_cellType(strcmp(responsivityVec_cellType,'-')), 0:0.04:1);
%     % hc(end+1) = 0;
%     % stairs(he, hc, 'Color', 'red', 'LineWidth',4);
% 
%     hold off
% 
%     xlabel('First Spike Latency (ms)');
%     ylabel('Unit Count');
% end


figure;

g = gramm('x',FSL_vec, 'color',groupVec, 'row',cellTypesVec);
g.stat_bin('edges',0:3:100, 'geom','overlaid_bar');
g.set_names('x','First Spike Latency (ms)', 'Color','', 'Row','');

g.draw;

end

