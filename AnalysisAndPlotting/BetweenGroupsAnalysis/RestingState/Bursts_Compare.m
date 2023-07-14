function [data_table_bursts, all_data] = Bursts_Compare(all_data, cell_types, plot_points)
% Plot and compare bursts between groups. Outputs a data table for stats,
% and updates all_data by adding some burst metrics.
%
% INPUTS
% - all_data
% - cell_types: cell array of cell types to plot. e.g., {'MSN','TAN'}
% - plot_points: 0 to only plot bars (+/- SEM), 1 to plot points on top.
%       Useful to show the distribution.

groupNames = fieldnames(all_data);

groupsVec = {};
cellTypesVec = {};
burstCount_vec = [];
propSpikesInBursts_vec = [];
proportion_time_in_bursts_vec = [];
mean_intra_burst_frequency_vec = [];
median_intra_burst_frequency_vec = [];
mean_spikes_per_burst_vec = [];
median_spikes_per_burst_vec = [];
bursts_per_second_vec = [];
mean_burst_duration_vec = [];
median_burst_duration_vec = [];

for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};

    recNames = fieldnames(all_data.(groupName));

    for recNum = 1:length(recNames)
        recName = recNames{recNum};

        cellIDs = fieldnames(all_data.(groupName).(recName));

        for cellID_num = 1:length(cellIDs)
            cellID = cellIDs{cellID_num};

            thisCellType = all_data.(groupName).(recName).(cellID).Cell_Type;
            isSingleUnit = all_data.(groupName).(recName).(cellID).IsSingleUnit;
            if any(strcmp(cell_types, thisCellType)) && isSingleUnit
                recLength = all_data.(groupName).(recName).(cellID).Recording_Duration;
                spikeTimes_samples = all_data.(groupName).(recName).(cellID).SpikeTimes_all;
                sr = all_data.(groupName).(recName).(cellID).Sampling_Frequency;
                spikeTimes_s = spikeTimes_samples/sr;

                burstsInfo = surpriseBurst(spikeTimes_s, recLength, 'min_length_of_burst',2, 'local_length',200, 'fac',1.5, 'surprise_cutoff',10);

                groupsVec{end+1,1} = groupName;
                cellTypesVec{end+1,1} = thisCellType;
                % metrics:
                burstCount_vec(end+1,1) = burstsInfo.num_bursts;
                propSpikesInBursts_vec(end+1,1) = burstsInfo.proportion_spikes_in_bursts;
                proportion_time_in_bursts_vec(end+1,1) = burstsInfo.proportion_time_in_bursts;
                mean_intra_burst_frequency_vec(end+1,1) = burstsInfo.mean_intra_burst_frequency;
                median_intra_burst_frequency_vec(end+1,1) = burstsInfo.median_intra_burst_frequency;
                mean_spikes_per_burst_vec(end+1,1) = burstsInfo.mean_spikes_per_burst;
                median_spikes_per_burst_vec(end+1,1) = burstsInfo.median_spikes_per_burst;
                bursts_per_second_vec(end+1,1) = burstsInfo.bursts_per_second;
                mean_burst_duration_vec(end+1,1) = burstsInfo.mean_burst_duration;
                median_burst_duration_vec(end+1,1) = burstsInfo.median_burst_duration;

                all_data.(groupName).(recName).(cellID).BurstsInfo = burstsInfo;
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
%         FRs_vec_sub = burstCount_vec(strcmp(groupsVec,groupName) & strcmp(cellTypesVec,cell_types{cell_type_ind}),1);
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
% num_bursts
figure;
g = gramm('x',groupsVec, 'y',burstCount_vec, 'color',groupsVec);
%g.facet_grid([],cellTypesVec, "scale","independent");
g.stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'width',1.2, 'setylim',true); % bar plot
g.set_names('x','', 'y','Number of Bursts', 'Color','', 'Column','');
g.no_legend;
g.draw();
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/num_bursts', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/num_bursts', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/num_bursts', 'png');
if plot_points % bar plot with points overlayed
    g.update('x',groupsVec, 'y',burstCount_vec, "color",groupsVec);
    g.geom_point();
    g.set_color_options('lightness',40);
    g.set_point_options("markers","^", "base_size",3);
    g.no_legend;
    g.draw;
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/num_bursts_withPoints', 'fig');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/num_bursts_withPoints', 'svg');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/num_bursts_withPoints', 'png');
end
figure;  % Probability density
g = gramm('x',burstCount_vec, 'color',groupsVec);
g.stat_density('npoints',1000);
g.set_names('x','Number of Bursts', 'y','Probability Density', 'color','');
g.draw;
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/num_bursts_PDF', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/num_bursts_PDF', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/num_bursts_PDF', 'png');

% proportion_spikes_in_bursts
figure;
g = gramm('x',groupsVec, 'y',propSpikesInBursts_vec, 'color',groupsVec);
%g.facet_grid([],cellTypesVec, "scale","independent");
g.stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'width',1.2, 'setylim',true);
g.set_names('x','', 'y','Proportion of Spikes in Bursts', 'Color','', 'Column','');
g.no_legend;
g.draw();
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/proportion_spikes_in_bursts', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/proportion_spikes_in_bursts', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/proportion_spikes_in_bursts', 'png');
if plot_points
    g.update('x',groupsVec, 'y',propSpikesInBursts_vec, "color",groupsVec);
    g.geom_point();
    g.set_color_options('lightness',40);
    g.set_point_options("markers","^", "base_size",3);
    g.no_legend;
    g.draw;
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/proportion_spikes_in_bursts_withPoints', 'fig');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/proportion_spikes_in_bursts_withPoints', 'svg');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/proportion_spikes_in_bursts_withPoints', 'png');
end
figure; % Probability density
g = gramm('x',propSpikesInBursts_vec, 'color',groupsVec);
g.stat_density('npoints',1000);
g.set_names('x','Proportion of Spikes in Bursts', 'y','Probability Density', 'color','');
g.draw;
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/proportion_spikes_in_bursts_PDF', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/proportion_spikes_in_bursts_PDF', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/proportion_spikes_in_bursts_PDF', 'png');

% proportion_time_in_bursts
figure;
g = gramm('x',groupsVec, 'y',proportion_time_in_bursts_vec, 'color',groupsVec);
%g.facet_grid([],cellTypesVec, "scale","independent");
g.stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'width',1.2, 'setylim',true);
g.set_names('x','', 'y','Proportion of Time in Bursts', 'Color','', 'Column','');
g.no_legend;
g.draw();
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/proportion_time_in_bursts', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/proportion_time_in_bursts', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/proportion_time_in_bursts', 'png');
if plot_points
    g.update('x',groupsVec, 'y',proportion_time_in_bursts_vec, "color",groupsVec);
    g.geom_point();
    g.set_color_options('lightness',40);
    g.set_point_options("markers","^", "base_size",3);
    g.no_legend;
    g.draw;
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/proportion_time_in_bursts_withPoints', 'fig');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/proportion_time_in_bursts_withPoints', 'svg');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/proportion_time_in_bursts_withPoints', 'png');
end
figure; % Probability density
g = gramm('x',proportion_time_in_bursts_vec, 'color',groupsVec);
g.stat_density('npoints',1000);
g.set_names('x','Proportion of Time in Bursts', 'y','Probability Density', 'color','');
g.draw;
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/proportion_time_in_bursts_PDF', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/proportion_time_in_bursts_PDF', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/proportion_time_in_bursts_PDF', 'png');

% mean_intra_burst_frequency
figure;
g = gramm('x',groupsVec, 'y',mean_intra_burst_frequency_vec, 'color',groupsVec);
%g.facet_grid([],cellTypesVec, "scale","independent");
g.stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'width',1.2, 'setylim',true);
g.set_names('x','', 'y','Mean Intra-Burst FR (Hz)', 'Color','', 'Column','');
g.no_legend;
g.draw();
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_intra_burst_frequency', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_intra_burst_frequency', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_intra_burst_frequency', 'png');
if plot_points
    g.update('x',groupsVec, 'y',mean_intra_burst_frequency_vec, "color",groupsVec);
    g.geom_point();
    g.set_color_options('lightness',40);
    g.set_point_options("markers","^", "base_size",3);
    g.no_legend;
    g.draw;
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_intra_burst_frequency_withPoints', 'fig');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_intra_burst_frequency_withPoints', 'svg');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_intra_burst_frequency_withPoints', 'png');
end
figure; % Probability density
g = gramm('x',mean_intra_burst_frequency_vec, 'color',groupsVec);
g.stat_density('npoints',1000);
g.set_names('x','Mean Intra-Burst FR (Hz)', 'y','Probability Density', 'color','');
g.draw;
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_intra_burst_frequency_PDF', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_intra_burst_frequency_PDF', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_intra_burst_frequency_PDF', 'png');

% median_intra_burst_frequency
figure;
g = gramm('x',groupsVec, 'y',median_intra_burst_frequency_vec, 'color',groupsVec);
%g.facet_grid([],cellTypesVec, "scale","independent");
g.stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'width',1.2, 'setylim',true);
g.set_names('x','', 'y','Median Intra-Burst FR (Hz)', 'Color','', 'Column','');
g.no_legend;
g.draw();
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_intra_burst_frequency', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_intra_burst_frequency', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_intra_burst_frequency', 'png');
if plot_points
    g.update('x',groupsVec, 'y',median_intra_burst_frequency_vec, "color",groupsVec);
    g.geom_point();
    g.set_color_options('lightness',40);
    g.set_point_options("markers","^", "base_size",3);
    g.no_legend;
    g.draw;
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_intra_burst_frequency_withPoints', 'fig');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_intra_burst_frequency_withPoints', 'svg');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_intra_burst_frequency_withPoints', 'png');
end
figure; % Probability density
g = gramm('x',median_intra_burst_frequency_vec, 'color',groupsVec);
g.stat_density('npoints',1000);
g.set_names('x','Median Intra-Burst FR (Hz)', 'y','Probability Density', 'color','');
g.draw;
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_intra_burst_frequency_PDF', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_intra_burst_frequency_PDF', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_intra_burst_frequency_PDF', 'png');

% mean_spikes_per_burst
figure;
g = gramm('x',groupsVec, 'y',mean_spikes_per_burst_vec, 'color',groupsVec);
%g.facet_grid([],cellTypesVec, "scale","independent");
g.stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'width',1.2, 'setylim',true);
g.set_names('x','', 'y','Mean Num. Spikes per Burst', 'Color','', 'Column','');
g.no_legend;
g.draw();
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_spikes_per_burst', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_spikes_per_burst', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_spikes_per_burst', 'png');
if plot_points
    g.update('x',groupsVec, 'y',mean_spikes_per_burst_vec, "color",groupsVec);
    g.geom_point();
    g.set_color_options('lightness',40);
    g.set_point_options("markers","^", "base_size",3);
    g.no_legend;
    g.draw;
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_spikes_per_burst_withPoints', 'fig');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_spikes_per_burst_withPoints', 'svg');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_spikes_per_burst_withPoints', 'png');
end
figure; % Probability density
g = gramm('x',mean_spikes_per_burst_vec, 'color',groupsVec);
g.stat_density('npoints',1000);
g.set_names('x','Mean Num. Spikes per Burst', 'y','Probability Density', 'color','');
g.draw;
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_spikes_per_burst_PDF', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_spikes_per_burst_PDF', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_spikes_per_burst_PDF', 'png');

% median_spikes_per_burst
figure;
g = gramm('x',groupsVec, 'y',median_spikes_per_burst_vec, 'color',groupsVec);
%g.facet_grid([],cellTypesVec, "scale","independent");
g.stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'width',1.2, 'setylim',true);
g.set_names('x','', 'y','Median Num. Spikes per Burst', 'Color','', 'Column','');
g.no_legend;
g.draw();
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_spikes_per_burst', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_spikes_per_burst', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_spikes_per_burst', 'png');
if plot_points
    g.update('x',groupsVec, 'y',median_spikes_per_burst_vec, "color",groupsVec);
    g.geom_point();
    g.set_color_options('lightness',40);
    g.set_point_options("markers","^", "base_size",3);
    g.no_legend;
    g.draw;
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_spikes_per_burst_withPoints', 'fig');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_spikes_per_burst_withPoints', 'svg');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_spikes_per_burst_withPoints', 'png');
end
figure; % Probability density
g = gramm('x',median_spikes_per_burst_vec, 'color',groupsVec);
g.stat_density('npoints',1000);
g.set_names('x','Median Num. Spikes per Burst', 'y','Probability Density', 'color','');
g.draw;
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_spikes_per_burst_PDF', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_spikes_per_burst_PDF', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_spikes_per_burst_PDF', 'png');

% bursts_per_second
figure;
g = gramm('x',groupsVec, 'y',bursts_per_second_vec, 'color',groupsVec);
%g.facet_grid([],cellTypesVec, "scale","independent");
g.stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'width',1.2, 'setylim',true);
g.set_names('x','', 'y','Bursts per Second', 'Color','', 'Column','');
g.no_legend;
g.draw();
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/bursts_per_second', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/bursts_per_second', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/bursts_per_second', 'png');
if plot_points
    g.update('x',groupsVec, 'y',bursts_per_second_vec, "color",groupsVec);
    g.geom_point();
    g.set_color_options('lightness',40);
    g.set_point_options("markers","^", "base_size",3);
    g.no_legend;
    g.draw;
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/bursts_per_second_withPoints', 'fig');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/bursts_per_second_withPoints', 'svg');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/bursts_per_second_withPoints', 'png');
end
figure; % Probability density
g = gramm('x',bursts_per_second_vec, 'color',groupsVec);
g.stat_density('npoints',1000);
g.set_names('x','Bursts per Second', 'y','Probability Density', 'color','');
g.draw;
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/bursts_per_second_PDF', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/bursts_per_second_PDF', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/bursts_per_second_PDF', 'png');

% mean_burst_duration
figure;
g = gramm('x',groupsVec, 'y',mean_burst_duration_vec, 'color',groupsVec);
%g.facet_grid([],cellTypesVec, "scale","independent");
g.stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'width',1.2, 'setylim',true);
g.set_names('x','', 'y','Mean Burst Duration (seconds)', 'Color','', 'Column','');
g.no_legend;
g.draw();
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_burst_duration', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_burst_duration', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_burst_duration', 'png');
if plot_points
    g.update('x',groupsVec, 'y',mean_burst_duration_vec, "color",groupsVec);
    g.geom_point();
    g.set_color_options('lightness',40);
    g.set_point_options("markers","^", "base_size",3);
    g.no_legend;
    g.draw;
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_burst_duration_withPoints', 'fig');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_burst_duration_withPoints', 'svg');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_burst_duration_withPoints', 'png');
end
figure; % Probability density
g = gramm('x',mean_burst_duration_vec, 'color',groupsVec);
g.stat_density('npoints',1000);
g.set_names('x','Mean Burst Duration (seconds)', 'y','Probability Density', 'color','');
g.draw;
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_burst_duration_PDF', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_burst_duration_PDF', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/mean_burst_duration_PDF', 'png');

% median_burst_duration
figure;
g = gramm('x',groupsVec, 'y',median_burst_duration_vec, 'color',groupsVec);
%g.facet_grid([],cellTypesVec, "scale","independent");
g.stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'width',1.2, 'setylim',true);
g.set_names('x','', 'y','Median Burst Duration (seconds)', 'Color','', 'Column','');
g.no_legend;
g.draw();
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_burst_duration', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_burst_duration', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_burst_duration', 'png');
if plot_points
    g.update('x',groupsVec, 'y',median_burst_duration_vec, "color",groupsVec);
    g.geom_point();
    g.set_color_options('lightness',40);
    g.set_point_options("markers","^", "base_size",3);
    g.no_legend;
    g.draw;
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_burst_duration_withPoints', 'fig');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_burst_duration_withPoints', 'svg');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_burst_duration_withPoints', 'png');
end
figure; % Probability density
g = gramm('x',median_burst_duration_vec, 'color',groupsVec);
g.stat_density('npoints',1000);
g.set_names('x','Median Burst Duration (seconds)', 'y','Probability Density', 'color','');
g.draw;
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_burst_duration_PDF', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_burst_duration_PDF', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Bursts/median_burst_duration_PDF', 'png');

%% Make data table to export for stats
data_table_bursts = table(groupsVec, cellTypesVec, burstCount_vec, propSpikesInBursts_vec, proportion_time_in_bursts_vec, mean_intra_burst_frequency_vec, median_intra_burst_frequency_vec, mean_spikes_per_burst_vec, median_spikes_per_burst_vec, bursts_per_second_vec, mean_burst_duration_vec, median_burst_duration_vec, ...
    'VariableNames',{'Group','CellType','num_bursts','proportion_spikes_in_bursts','proportion_time_in_bursts','mean_intra_burst_frequency','median_intra_burst_frequency','mean_spikes_per_burst','median_spikes_per_burst','bursts_per_second','mean_burst_duration','median_burst_duration'});

end

