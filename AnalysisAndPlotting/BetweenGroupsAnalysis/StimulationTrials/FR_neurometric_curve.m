function stats_out = FR_neurometric_curve(all_data,trialTagsLabels,userGroupNames,groupsPlottingOrder,FR_type,excludeMUA)
% Plot a neurometric curve showing the mean response of each cell type to
% varying levels of stimulation (e.g., whisker stim intensities) between
% groups. Also runs stats and returns the results as a struct.
%
% INPUTS
%   userGroupNames: {'group1','group2',...} to rename groups in the same
%       order that they appear in the all_data struct.
%   groupsPlottingOrder: {'group1','group2',...}, with the names given in
%       userGroupNames, to determine the order of plotting in Gramm.
%   FR_type: 'binned', 'inst', 'peak', or 'fano' (not really a type of FR
%       but ya know, it's useful)
%   excludeMUA: 0 or 1


groupNames = fieldnames(all_data);

FRsVec = []; % vector with row for each unit
groupsVec = {};
trialTagsVec = {};
cellTypesVec = {};
cellTypesVec2 = {};
responsivityVec = {};
layerVec = {};

for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};
    userGroupName = userGroupNames{groupNum};
    mouseNames = fieldnames(all_data.(groupName));

    for mouseNum = 1:length(mouseNames)
        mouseName = mouseNames{mouseNum};
        cellIDs = fieldnames(all_data.(groupName).(mouseName));

        for cellID_num = 1:length(cellIDs)
            cellID = cellIDs{cellID_num};

            IsSingleUnit = all_data.(groupName).(mouseName).(cellID).IsSingleUnit;
            ISI_violations_percent = all_data.(groupName).(mouseName).(cellID).ISI_violations_percent;
            
            %if (excludeMUA && IsSingleUnit) || ~excludeMUA
            if (excludeMUA && (ISI_violations_percent <= 1)) || ~excludeMUA
                for trialTagInd = 1:length(trialTagsLabels)
                    
                    if strcmp(FR_type,'binned')
                        FRsVec(end+1,1) = all_data.(groupName).(mouseName).(cellID).MeanFR_stim(trialTagInd,1);
                    elseif strcmp(FR_type,'inst')
                        FRsVec(end+1,1) = all_data.(groupName).(mouseName).(cellID).MeanFR_inst_stim(trialTagInd,1);
                    elseif strcmp(FR_type,'peak')
                        FRsVec(end+1,1) = all_data.(groupName).(mouseName).(cellID).PeakEvokedFR(trialTagInd,1);
                    else % FR_type is 'fano'
                        FRsVec(end+1,1) = all_data.(groupName).(mouseName).(cellID).FanoFactor_stim(trialTagInd,1);
                    end
    
                    cell_type = all_data.(groupName).(mouseName).(cellID).Cell_Type;

                    trialTagsVec{end+1,1} = trialTagsLabels{trialTagInd};
                    cellTypesVec{end+1,1} = cell_type;
                    groupsVec{end+1,1} = userGroupName;
                    layerVec{end+1,1} = all_data.(groupName).(mouseName).(cellID).LaminarLabel;
    
                    responsivityNum = all_data.(groupName).(mouseName).(cellID).StimResponsivity;
                    if responsivityNum == 1
                        responsivityLabel = '+';
                        responsivityVec{end+1,1} = responsivityLabel;
                    elseif responsivityNum == 0
                        responsivityLabel = 'nr';
                        responsivityVec{end+1,1} = responsivityLabel;
                    else
                        responsivityLabel = '-';
                        responsivityVec{end+1,1} = responsivityLabel;
                    end

                    cellTypesVec2{end+1,1} = strcat(responsivityLabel,cell_type);
                    
                end
            end

        end
    end
end

cell_types = unique(cellTypesVec);

figure();

%% Plotting response curves
g = gramm('x',trialTagsVec, 'y',FRsVec, 'color',groupsVec, 'column',cellTypesVec2, 'row',layerVec);
g.stat_summary('type','sem', 'geom',{'line','black_errorbar','point'}, 'setylim',true);
if strcmp(FR_type,'peak')
    g.set_names('x','', 'y','Peak Evoked FR (Hz)', 'Color','', 'Column','', 'Row','');
elseif strcmp(FR_type,'fano')
    g.set_names('x','', 'y','Fano Factor', 'Color','', 'Column','', 'Row','');
else
    g.set_names('x','', 'y','Firing Rate (Hz)', 'Color','', 'Column','', 'Row','');
end
g.set_order_options('x',trialTagsLabels, "color",groupsPlottingOrder, 'row',{'SG','L4','IG'});
g.draw();

%% Generate table of sample sizes
% for ii = 1:length(unique(cellTypesVec)) % rows = cell types
%     for jj = 1:length(unique(responsivityVec)) % columns = responsivity
%         for k = 1:length(unique(groupsVec)) % subcolumns = groups
%             n_groupsVec = 
%         end
%     end
% end
% 
% g2 = gramm('x',groupsVec, )

%% Stats: repeated measures ANOVA and Bonferroni test for multiple
% comparisons
responsivityLabels = unique(responsivityVec);

group1_inds = find(strcmp(groupsVec,userGroupNames{1}));
group2_inds = find(strcmp(groupsVec,userGroupNames{2}));
for uniqueCellTypeInd = 1:length(cell_types)
    withinDsgn = table((0:length(trialTagsLabels)-1)','VariableNames',{'Intensity'});
    cellTypeInds = find(strcmp(cellTypesVec,cell_types{uniqueCellTypeInd}));

    for responsivityLabelInd = 1:length(responsivityLabels)
        responsivityLabel = responsivityLabels{responsivityLabelInd};
        thisResponsivityLabelInds = find(strcmp(responsivityVec,responsivityLabel));

        % '+' and '-' aren't valid fieldnames for the output struct, so
        % replace them with words
        if strcmp(responsivityLabel,'+')
            responsivityLabel = 'Positive';
        elseif strcmp(responsivityLabel,'-')
            responsivityLabel = 'Negative';
        end

        for trialTagLabelInd = 1:length(trialTagsLabels)
            trialTagLabel = trialTagsLabels{trialTagLabelInd};
            trialTagInds = find(strcmp(trialTagsVec,trialTagLabel));
    
            cellType_responsivityLabel_inds = intersect(cellTypeInds,thisResponsivityLabelInds);
            cellType_responsivityLabel_trialTag_inds = intersect(cellType_responsivityLabel_inds,trialTagInds);
            cellType_responsivityLabel_trialTag_group1_inds = intersect(cellType_responsivityLabel_trialTag_inds,group1_inds);
            cellType_responsivityLabel_trialTag_group2_inds = intersect(cellType_responsivityLabel_trialTag_inds,group2_inds);
    
            FRs_group1_trialTag_responsivityLabel_cellType = FRsVec(cellType_responsivityLabel_trialTag_group1_inds);
            FRs_group2_trialTag_responsivityLabel_cellType = FRsVec(cellType_responsivityLabel_trialTag_group2_inds);
    
            [h,p,ci,stats] = ttest2(FRs_group1_trialTag_responsivityLabel_cellType, FRs_group2_trialTag_responsivityLabel_cellType, 'Vartype','unequal');
    
            n_table = table(length(FRs_group1_trialTag_responsivityLabel_cellType),length(FRs_group2_trialTag_responsivityLabel_cellType), 'VariableNames',userGroupNames);

            stats_out.(cell_types{uniqueCellTypeInd}).(responsivityLabel).(trialTagLabel).n = n_table;
            stats_out.(cell_types{uniqueCellTypeInd}).(responsivityLabel).(trialTagLabel).h = h;
            stats_out.(cell_types{uniqueCellTypeInd}).(responsivityLabel).(trialTagLabel).p = p;
            stats_out.(cell_types{uniqueCellTypeInd}).(responsivityLabel).(trialTagLabel).ci = ci;
            stats_out.(cell_types{uniqueCellTypeInd}).(responsivityLabel).(trialTagLabel).stats = stats;
        end
    end
end

end

