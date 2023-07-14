function data_table_FR = FR_compare(all_data, cell_types, binSize, plot_points)
% Plot and compare ISI coefficient of variation (CV) between groups.
%
% INPUTS
% - all_data
% - cell_types: cell array of cell types to plot. e.g., {'MSN','TAN'}
% - binSize: size of bins (seconds) to calculate FRs. The maximum of those
%            binned FRs will be the final saved FR of the unit. If 0, the
%            whole-recording FR is used instead.
% - plot_points: 0 to only plot bars (+/- SEM), 1 to plot points on top.
%       Useful to show the distribution.

groupNames = fieldnames(all_data);

groupsVec = {};
cellTypesVec = {};
FRs_vec = [];

for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};

    mouseNames = fieldnames(all_data.(groupName));

    for mouseNum = 1:length(mouseNames)
        mouseName = mouseNames{mouseNum};

        cellIDs = fieldnames(all_data.(groupName).(mouseName));

        for cellID_num = 1:length(cellIDs)
            cellID = cellIDs{cellID_num};

            thisCellType = all_data.(groupName).(mouseName).(cellID).Cell_Type;
            isSingleUnit = all_data.(groupName).(mouseName).(cellID).IsSingleUnit;
            if any(strcmp(cell_types, thisCellType)) && isSingleUnit
                groupsVec{end+1,1} = groupName;
                cellTypesVec{end+1,1} = thisCellType;

                if binSize == 0
                    FRs_vec(end+1,1) = all_data.(groupName).(mouseName).(cellID).MeanFR_total;
                else
                    % calculate FR as max of 200-second bin FRs
                    spikeTimes = all_data.(groupName).(mouseName).(cellID).SpikeTimes_all / all_data.(groupName).(mouseName).(cellID).Sampling_Frequency;
                    RecDuration = all_data.(groupName).(mouseName).(cellID).Recording_Duration;
                    intervalBounds = 0:binSize:RecDuration;
                    binned_FRs_vec = [];
                    for ii = 1:length(intervalBounds)-1
                        n_spikes = length(spikeTimes((spikeTimes >= intervalBounds(ii))&(spikeTimes <= intervalBounds(ii+1))));
                        binned_FRs_vec(end+1,1) = n_spikes / binSize;
                    end
                    FRs_vec(end+1,1) = max(binned_FRs_vec);
                end
            end
        end
    end
end

%% Remove outliers
% FRs_vec_new = [];
% groupsVec_new = {};
% cellTypesVec_new = {};
% 
% for groupNum = 1:length(groupNames)
%     groupName = groupNames{groupNum};
%     for cell_type_ind = 1:length(cell_types)
%         cell_type_name = cell_types{cell_type_ind};
% 
%         FRs_vec_sub = FRs_vec(strcmp(groupsVec,groupName) & strcmp(cellTypesVec,cell_types{cell_type_ind}),1);
%         FRs_vec_sub = rmoutliers(FRs_vec_sub);
% 
%         groupsVec_sub = groupsVec(strcmp(groupsVec,groupName),1);
%         groupsVec_sub = groupsVec_sub(1:length(FRs_vec_sub),1);
% 
%         cellTypesVec_sub = cellTypesVec(strcmp(cellTypesVec,cell_type_name),1);
%         cellTypesVec_sub = cellTypesVec_sub(1:length(FRs_vec_sub),1);
% 
%         FRs_vec_new = [FRs_vec_new; FRs_vec_sub];
%         groupsVec_new = cat(1, groupsVec_new, groupsVec_sub);
%         cellTypesVec_new = cat(1, cellTypesVec_new, cellTypesVec_sub);
%     end
% end
% FRs_vec = FRs_vec_new;
% groupsVec = groupsVec_new;
% cellTypesVec = cellTypesVec_new;

%% Plotting
figure;
g = gramm('x',groupsVec, 'y',FRs_vec, 'color',groupsVec);
%g.facet_grid([],cellTypesVec, "scale","independent");
g.stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'width',1.2, 'setylim',true);
g.set_names('x','', 'y','Firing Rate (Hz)', 'Color','', 'Column','');
g.no_legend;
g.draw();
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/FiringRate/FR_barplot', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/FiringRate/FR_barplot', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/FiringRate/FR_barplot', 'png');
if plot_points
    g.update('x',groupsVec, 'y',FRs_vec, "color",groupsVec);
    g.geom_point();
    g.set_color_options('lightness',40);
    g.set_point_options("markers","^", "base_size",3);
    g.no_legend;
    g.draw;
end
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/FiringRate/FR_barplot_withPoints', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/FiringRate/FR_barplot_withPoints', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/FiringRate/FR_barplot_withPoints', 'png');

%% Make data table to export for stats
data_table_FR = table(groupsVec, cellTypesVec, FRs_vec, 'VariableNames',{'Group','CellType','FR'});

end

