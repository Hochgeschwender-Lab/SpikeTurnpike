function ISI_CV_by_group_histogram(all_data)
%AUROC_HISTOGRAM Summary of this function goes here
%   Detailed explanation goes here

groupNames = fieldnames(all_data);

CV_vec = [];
cellTypesVec = {};
responsivityVec = {};
groupVec = {};

for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};
    mouseNames = fieldnames(all_data.(groupName));

    for mouseNum = 1:length(mouseNames)
        mouseName = mouseNames{mouseNum};
        cellIDs = fieldnames(all_data.(groupName).(mouseName));

        for cellID_num = 1:length(cellIDs)
            cellID = cellIDs{cellID_num};

            CV_vec(end+1,1) = all_data.(groupName).(mouseName).(cellID).ISI_baseline_CV;
            cellTypesVec{end+1,1} = all_data.(groupName).(mouseName).(cellID).Cell_Type;
            groupVec{end+1,1} = groupName;

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

% figure;
% T = tiledlayout('flow');
% for cellTypeNum = 1:length(cell_types)
%     CV_vec_cellType = CV_vec(strcmp(cellTypesVec, cell_types{cellTypeNum}));
%     responsivityVec_cellType = responsivityVec(strcmp(cellTypesVec, cell_types{cellTypeNum}));
% 
%     t = nexttile(T);
%     title(t, cell_types{cellTypeNum});
%     hold on
%     h1 = histogram(CV_vec_cellType(strcmp(responsivityVec_cellType,'+')), 'BinWidth',0.1);
%     h2 = histogram(CV_vec_cellType(strcmp(responsivityVec_cellType,'NR')), 'BinWidth',0.1);
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
%     xlabel('Coefficient of Variation');
%     ylabel('Unit Count');
% end


figure;

g = gramm('x',CV_vec, 'color',groupVec, 'row',cellTypesVec);
g.stat_bin('edges',0:0.1:3, 'geom','overlaid_bar');
g.set_names('x','Coefficient of Variation', 'Color','', 'Row','');

g.draw;

end

