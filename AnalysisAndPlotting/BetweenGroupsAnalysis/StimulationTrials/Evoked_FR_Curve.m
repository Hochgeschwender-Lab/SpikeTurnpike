function stats_out = Evoked_FR_Curve(all_data,trialTagsLabels,userGroupNames,groupsPlottingOrder,FR_type,minNspikes,excludeMUA)
%EVOKED_FR_CURVE Summary of this function goes here
%   
%   Summary
%       This function plots the evoked firing rate (FR) curve and performs statistical analyses for given neural data. Depending on the trial tags provided, it can either generate bar plots for invariant stimulation or response curves for varying stimulation. Statistical tests, including two-sample t-tests and repeated measures ANOVA, are applied to the data.
%   
%   Input Parameters
%       all_data: A structured array containing the detailed information for different groups, mice, and cell IDs. The structure includes fields like MeanFR_stim, MeanFR_baseline, Cell_Type, and StimResponsivity.
%       trialTagsLabels: A cell array containing the labels for the trial tags. It defines the stimulation conditions to be plotted and analyzed.
%       userGroupNames: A cell array containing the names of user-defined groups to be used in the plot.
%       groupsPlottingOrder: A cell array defining the order in which the groups should be plotted.
%   Output
%       stats_out: A structured array containing the statistical analysis results. It includes the p-values, confidence intervals, t-statistic, and degrees of freedom for the conducted tests.

groupNames = fieldnames(all_data);

FRsVec = []; % vector with row for each unit
groupsVec = {};
trialTagsVec = {}; % numeric trial tag, not label
cellTypesVec = {};
responsivityVec = {};
cellTypesVec2 = {};

for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};
    userGroupName = userGroupNames{groupNum};
    mouseNames = fieldnames(all_data.(groupName));

    for mouseNum = 1:length(mouseNames)
        mouseName = mouseNames{mouseNum};
        cellIDs = fieldnames(all_data.(groupName).(mouseName));

        for cellID_num = 1:length(cellIDs)
            cellID = cellIDs{cellID_num};
            ISI_violations_percent = all_data.(groupName).(mouseName).(cellID).ISI_violations_percent;
            nspikes = length(all_data.(groupName).(mouseName).(cellID).SpikeTimes_all);

            if ((excludeMUA && (ISI_violations_percent <= 1.5)) || ~excludeMUA) && (nspikes >= minNspikes)
                if length(trialTagsLabels) > 1
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
                        trialTagsVec{end+1,1} = trialTagsLabels{trialTagInd};
                        cell_type = all_data.(groupName).(mouseName).(cellID).Cell_Type;
                        cellTypesVec{end+1,1} = cell_type;
                        groupsVec{end+1,1} = userGroupName;
    
                        responsivityNum = all_data.(groupName).(mouseName).(cellID).StimResponsivity;
                        if responsivityNum == 1
                            responsivityLabel = '+';
                            responsivityVec{end+1,1} = 'Positive';
                        elseif responsivityNum == 0
                            responsivityLabel = 'nr';
                            responsivityVec{end+1,1} = 'NR';
                        else
                            responsivityLabel = '-';
                            responsivityVec{end+1,1} = 'Negative';
                        end

                        cellTypesVec2{end+1,1} = strcat(responsivityLabel,cell_type);
                        
                    end
                else
                    FRsVec(end+1,1) = all_data.(groupName).(mouseName).(cellID).MeanFR_baseline;
                    trialTagsVec{end+1,1} = 'Baseline';
                    cell_type = all_data.(groupName).(mouseName).(cellID).Cell_Type;
                    cellTypesVec{end+1,1} = cell_type;
                    groupsVec{end+1,1} = userGroupName;
                    responsivityNum = all_data.(groupName).(mouseName).(cellID).StimResponsivity;
                    if responsivityNum == 1
                        responsivityLabel = '+';
                        responsivityVec{end+1,1} = 'Positive';
                    elseif responsivityNum == 0
                        responsivityLabel = 'nr';
                        responsivityVec{end+1,1} = 'NR';
                    else
                        responsivityLabel = '-';
                        responsivityVec{end+1,1} = 'Negative';
                    end
                    cellTypesVec2{end+1,1} = strcat(responsivityLabel,cell_type);
                    
                    FRsVec(end+1,1) = all_data.(groupName).(mouseName).(cellID).MeanFR_stim;
                    trialTagsVec{end+1,1} = 'Evoked';
                    cell_type = all_data.(groupName).(mouseName).(cellID).Cell_Type;
                    cellTypesVec{end+1,1} = cell_type;
                    groupsVec{end+1,1} = userGroupName;
                    responsivityVec{end+1,1} = all_data.(groupName).(mouseName).(cellID).StimResponsivity;
                    responsivityNum = all_data.(groupName).(mouseName).(cellID).StimResponsivity;
                    if responsivityNum == 1
                        responsivityLabel = '+';
                        responsivityVec{end+1,1} = 'Positive';
                    elseif responsivityNum == 0
                        responsivityLabel = 'nr';
                        responsivityVec{end+1,1} = 'NR';
                    else
                        responsivityLabel = '-';
                        responsivityVec{end+1,1} = 'Negative';
                    end
                    cellTypesVec2{end+1,1} = strcat(responsivityLabel,cell_type);
                end
            end
        end
    end
end

cell_types = unique(cellTypesVec);

figure();

if length(trialTagsLabels) == 1 % bar plot + stats for invariant stimulation
    g = gramm('x',groupsVec, 'y',FRsVec, 'color',groupsVec);
    g.facet_grid(cellTypesVec2, trialTagsVec, 'scale','free_y');
    g.stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'setylim',true);
    g.set_names('x','', 'y','Firing Rate (Hz)', 'Color','', 'Column','', 'Row','');
    g.no_legend;
    g.set_order_options('x',groupsPlottingOrder, 'color',groupsPlottingOrder);
    g.draw();

    g.update('x',groupsVec, 'y',FRsVec, 'color',groupsVec);
    g.geom_jitter();
    g.set_color_options('lightness',40);
    g.set_point_options("markers","^", "base_size",3);
    g.no_legend;
    g.draw;
    
    for uniqueCellTypeInd = 1:length(cell_types)
        cellTypeInds = find(strcmp(cellTypesVec,cell_types{uniqueCellTypeInd}));
        % Stats for baseline
        baseline_inds = find(strcmp(trialTagsVec,'Baseline'));
        group1_inds = find(strcmp(groupsVec,userGroupNames{1}));
        group1_baseline_inds = intersect(baseline_inds,group1_inds);
        group1_baseline_cellType_inds = intersect(group1_baseline_inds,cellTypeInds);
        group2_inds = find(strcmp(groupsVec,userGroupNames{2}));
        group2_baseline_inds = intersect(baseline_inds,group2_inds);
        group2_baseline_cellType_inds = intersect(group2_baseline_inds,cellTypeInds);
        FRs_group1_baseline_cellType = FRsVec(group1_baseline_cellType_inds);
        FRs_group2_baseline_cellType = FRsVec(group2_baseline_cellType_inds);
        % Two-sample t-test for unequal variances
        [h_baseline,p_baseline,ci_baseline,stats_baseline] = ttest2(FRs_group1_baseline_cellType, FRs_group2_baseline_cellType, 'Vartype','unequal');
        stats_out.Baseline.(cell_types{uniqueCellTypeInd}).h = h_baseline;
        stats_out.Baseline.(cell_types{uniqueCellTypeInd}).p = p_baseline;
        stats_out.Baseline.(cell_types{uniqueCellTypeInd}).ci = ci_baseline;
        stats_out.Baseline.(cell_types{uniqueCellTypeInd}).stats = stats_baseline;
    
        % Stats for stimulation window
        evoked_inds = find(strcmp(trialTagsVec,'Evoked'));
        group1_evoked_inds = intersect(evoked_inds,group1_inds);
        group1_evoked_cellType_inds = intersect(group1_evoked_inds,cellTypeInds);
        group2_evoked_inds = intersect(evoked_inds,group2_inds);
        group2_evoked_cellType_inds = intersect(group2_evoked_inds,cellTypeInds);
        FRs_group1_evoked_cellType = FRsVec(group1_evoked_cellType_inds);
        FRs_group2_evoked_cellType = FRsVec(group2_evoked_cellType_inds);
        [h_evoked,p_evoked,ci_evoked,stats_evoked] = ttest2(FRs_group1_evoked_cellType, FRs_group2_evoked_cellType, 'Vartype','unequal');
        stats_out.Evoked.(cell_types{uniqueCellTypeInd}).h = h_evoked;
        stats_out.Evoked.(cell_types{uniqueCellTypeInd}).p = p_evoked;
        stats_out.Evoked.(cell_types{uniqueCellTypeInd}).ci = ci_evoked;
        stats_out.Evoked.(cell_types{uniqueCellTypeInd}).stats = stats_evoked;
    end

else % response curves + stats for varying stimulation
    g = gramm('x',trialTagsVec, 'y',FRsVec, 'color',groupsVec, 'column',cellTypesVec, 'row',responsivityVec);
    g.stat_summary('type','sem', 'geom',{'line','black_errorbar','point'}, 'setylim',true);
    g.set_names('x','', 'y','Firing Rate (Hz)', 'Color','', 'Column','')
    g.set_order_options('x',trialTagsLabels, "color",groupsPlottingOrder);
    g.draw();

    % Stats: repeated measures ANOVA and Bonferroni test for multiple
    % comparisons
    group1_inds = find(strcmp(groupsVec,userGroupNames{1}));
    group2_inds = find(strcmp(groupsVec,userGroupNames{2}));
    for uniqueCellTypeInd = 1:length(cell_types)
        withinDsgn = table((0:length(trialTagsLabels)-1)','VariableNames',{'Intensity'});
        cellTypeInds = find(strcmp(cellTypesVec,cell_types{uniqueCellTypeInd}));

        for trialTagLabelInd = 1:length(trialTagsLabels)
            trialTagLabel = trialTagsLabels{trialTagLabelInd};
            trialTagInds = find(strcmp(trialTagsVec,trialTagLabel));

            cellType_trialTag_inds = intersect(cellTypeInds,trialTagInds);
            cellType_trialTag_group1_inds = intersect(cellType_trialTag_inds,group1_inds);
            cellType_trialTag_group2_inds = intersect(cellType_trialTag_inds,group2_inds);

            FRs_group1_trialTag_cellType = FRsVec(cellType_trialTag_group1_inds);
            FRs_group2_trialTag_cellType = FRsVec(cellType_trialTag_group2_inds);

            [h,p,ci,stats] = ttest2(FRs_group1_trialTag_cellType, FRs_group2_trialTag_cellType, 'Vartype','unequal');

            stats_out.(cell_types{uniqueCellTypeInd}).(trialTagLabel).h = h;
            stats_out.(cell_types{uniqueCellTypeInd}).(trialTagLabel).p = p;
            stats_out.(cell_types{uniqueCellTypeInd}).(trialTagLabel).ci = ci;
            stats_out.(cell_types{uniqueCellTypeInd}).(trialTagLabel).stats = stats;
        end
end

end

