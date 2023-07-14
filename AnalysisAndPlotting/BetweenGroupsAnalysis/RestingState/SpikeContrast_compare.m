function [data_table_syn, s_curves, s_bins] = SpikeContrast_compare(all_data, cell_types, plot_points)
% Per cell type, calculates synchrony via the Spike-Contrast method to
% produce two plots: a bar plot of synchrony (S) between groups as well as
% the mean synchrony curve per group.
%
% INPUTS
% - all_data
% - cell_types: cell array of cell types to plot. e.g., {'MSN','TAN'}.
%       Unlike in other functions, the cell types are pooled for analysis.
% - plot_points: 0 to only plot bars (+/- SEM), 1 to plot points on top.
%       Useful to show the distribution.

groupNames = fieldnames(all_data);

% vectors to be populated: all should have length numRecordings
S_vec = [];
groupsVec = {};

s_curves = []; % synchrony curves, [numRecordings x nBins]

for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};

    recNames = fieldnames(all_data.(groupName));

    for recNum = 1:length(recNames)
        recName = recNames{recNum};

        fprintf('  Working on recording %s\n', recName);

        cellIDs = fieldnames(all_data.(groupName).(recName));

        % construct allSpikeTrains array
        allSpikeTrains = []; % [maxNumOfSpikes x units], must be NaN-padded
        for cellID_num = 1:length(cellIDs)
            cellID = cellIDs{cellID_num};

            thisCellType = all_data.(groupName).(recName).(cellID).Cell_Type;
            isSingleUnit = all_data.(groupName).(recName).(cellID).IsSingleUnit;

            if any(strcmp(cell_types, thisCellType)) && isSingleUnit %&& (all_data.(groupName).(recName).(cellID).MeanFR_total < 5)
                spikeTimes = all_data.(groupName).(recName).(cellID).SpikeTimes_all'; % transposed so it's a row vector
                spikeTimes = spikeTimes / 30000; % convert from samples to seconds
                n_spikes = length(spikeTimes);

                if isempty(allSpikeTrains) || (n_spikes == length(allSpikeTrains))
                    allSpikeTrains(:,end+1) = spikeTimes;
                    recLength = all_data.(groupName).(recName).(cellID).Recording_Duration; % getting recording duration here so it's only done once
                elseif n_spikes > length(allSpikeTrains)
                    allSpikeTrains(end+1:n_spikes,:) = NaN; % NaN-pad allSpikeTrains to the length of this spike train
                    allSpikeTrains(:,end+1) = spikeTimes;
                else % n_spikes < length(allSpikeTrains)
                    spikeTimes(end+1:length(allSpikeTrains)) = NaN; % NaN-pad this spike train to the length of allSpikeTrains
                    allSpikeTrains(:,end+1) = spikeTimes;
                end
            end
        end

        % Spike-Contrast implementation
        [S,PREF] = SpikeContrast(allSpikeTrains, recLength, 0.001, 2); % S is the synchrony index, 0 <= S <= 1
        s_curve = PREF.s;
        s_bins = PREF.bins;

        S_vec(end+1,1) = S;
        s_curves(end+1,:) = s_curve;
        groupsVec{end+1,1} = groupName;

        % % plot synchrony curve per recording (sanity check)
        % figure;
        % plot(s_bins, s_curve);
        % title(recName);

    end
end

%% Bar plot
figure;

g = gramm('x',groupsVec, 'y',S_vec, 'color',groupsVec);
g.stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'width',1.2, 'setylim',true);
g.set_names('x','', 'y','Synchrony Index', 'Color','');
g.no_legend;
g.draw();
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Spike-Contrast/SpikeContrast_BarPlot', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Spike-Contrast/SpikeContrast_BarPlot', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Spike-Contrast/SpikeContrast_BarPlot', 'png');
if plot_points
    g.update('x',groupsVec, 'y',S_vec, "color",groupsVec);
    g.geom_point();
    g.set_color_options('lightness',40);
    g.set_point_options("markers","^", "base_size",3);
    g.no_legend;
    g.draw;
end
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Spike-Contrast/SpikeContrast_BarPlot_withPoints', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Spike-Contrast/SpikeContrast_BarPlot_withPoints', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Spike-Contrast/SpikeContrast_BarPlot_withPoints', 'png');

%% Mean synchrony curves
figure;
g = gramm('x',s_bins, 'y',s_curves, 'color',groupsVec);
g.stat_summary('type','sem', 'setylim',true);
g.set_names('x','Time bin (s)', 'y','Synchrony', 'Color','');
g.draw;
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Spike-Contrast/SynchronyCurves_mean', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Spike-Contrast/SynchronyCurves_mean', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Spike-Contrast/SynchronyCurves_mean', 'png');

%% Individual synchrony curves
figure;
g = gramm('x',s_bins, 'y',s_curves, 'color',groupsVec);
g.geom_line();
g.set_names('x','Time bin (s)', 'y','Synchrony', 'Color','');
g.draw;
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Spike-Contrast/SynchronyCurves_all', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Spike-Contrast/SynchronyCurves_all', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/Spike-Contrast/SynchronyCurves_all', 'png');

%% Save data table
data_table_syn = table(groupsVec, S_vec, 'VariableNames',{'Group','Synchrony'});

end