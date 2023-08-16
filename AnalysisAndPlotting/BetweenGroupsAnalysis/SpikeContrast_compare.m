function [data_table_syn, s_curves, s_bins] = SpikeContrast_compare(all_data, trialTagsLabels, binSizes, minNspikes, plot_points)
% Per cell type, calculates synchrony via the Spike-Contrast method to
% produce two plots: a bar plot of synchrony (S) between groups as well as
% the mean synchrony curve per group.
%
% INPUTS
% - all_data
% - trialTagsLabels: e.g., {'Catch', 'Low', 'Mid', 'Max'}
% - binSizes: [min max] bin sizes (seconds), i.e. timescales, for synchrony
%       calculation. WARNING: max bin size should not exceed half of your
%       recording or trial length, or the results will not make sense.
% - minNspikes: minimum number of times a unit must spike _within a trial
%       condition_ to be used for the SpikeContrast _of that trial type_.
% - plot_points: 0 to only plot bars (+/- SEM), 1 to plot points on top.
%       Useful to show the distribution.

groupNames = fieldnames(all_data);

% vectors (and one array) to be populated: all should have length numRecordings
S_vec = [];
s_curves = []; % synchrony curves, [numRecordings x nBins]
groupsVec = {};
cellTypesVec = {};
trialTypesVec = {};
%cellTypesVec2 = {};
responsivityVec = {};

for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};

    recNames = fieldnames(all_data.(groupName));

    for recNum = 1:length(recNames)
        recName = recNames{recNum};

        fprintf('  Working on recording %s\n', recName);

        cellIDs = fieldnames(all_data.(groupName).(recName));

        % construct allSpikeTrains array
        allSpikeTrains = []; % [maxNumOfSpikes x units], must be NaN-padded
        cellTypesVec_spiketrains = {}; % keep track of cell types and trial types for indexing allSpikeTrains
        %cellTypesVec2_spiketrains = {};
        ResponsivityVec_spiketrains = {};
        trialTypesVec_spiketrains = {};
        for cellID_num = 1:length(cellIDs)
            cellID = cellIDs{cellID_num};

            cell_type = all_data.(groupName).(recName).(cellID).Cell_Type;
            StimResp = all_data.(groupName).(recName).(cellID).StimResponsivity;
            if StimResp == 1
                StimResp = '+';
            elseif StimResp == 0
                StimResp = 'NR';
            elseif StimResp == -1
                StimResp = '-';
            end
            %cell_type2 = strcat(StimResp,cell_type);
            
            % necessary because our numeric label for '8 Hz LED' was 8
            % instead of 1 for some reason
            if groupNum==1 && recNum==1 && cellID_num==1
                trialTagsLabels_numeric = unique(all_data.(groupName).(recName).(cellID).Stim_Intensity);
            end

            %isSingleUnit = all_data.(groupName).(recName).(cellID).IsSingleUnit;
            ISI_violations_percent = all_data.(groupName).(recName).(cellID).ISI_violations_percent;

            if (ISI_violations_percent <= 1.5) %&& (all_data.(groupName).(recName).(cellID).MeanFR_total < 5)
                for trialType_ind = 1:length(trialTagsLabels)
                    trial_inds = all_data.(groupName).(recName).(cellID).Stim_Intensity == trialTagsLabels_numeric(trialType_ind);
                    spikeTimes = cell2mat(all_data.(groupName).(recName).(cellID).SpikeTimes_stim(1,trial_inds))'; % transposed so it's a row vector
                    spikeTimes = spikeTimes / 30000; % convert from samples to seconds
    
                    n_spikes = length(spikeTimes);
                    
                    if n_spikes >= minNspikes
                        % this complicated-looking block just handles NaN
                        % padding while constructing the allSpikeTrains array
                        if isempty(allSpikeTrains) || (n_spikes == size(allSpikeTrains,1))
                            allSpikeTrains(:,end+1) = spikeTimes;
                            recLength = all_data.(groupName).(recName).(cellID).Recording_Duration; % getting recording duration here so it's only done once
                        elseif n_spikes > size(allSpikeTrains,1)
                            allSpikeTrains(end+1:n_spikes,:) = NaN; % NaN-pad allSpikeTrains to the length of this spike train
                            allSpikeTrains(:,end+1) = spikeTimes;
                        else % n_spikes < size(allSpikeTrains,1)
                            spikeTimes(end+1:size(allSpikeTrains,1)) = NaN; % NaN-pad this spike train to the length of allSpikeTrains
                            allSpikeTrains(:,end+1) = spikeTimes;
                        end
    
                        cellTypesVec_spiketrains{1,end+1} = cell_type;
                        %cellTypesVec2_spiketrains{end+1,1} = cell_type2;
                        ResponsivityVec_spiketrains{1,end+1} = StimResp;
                        trialTypesVec_spiketrains{1,end+1} = trialTagsLabels{trialType_ind};
                    end
                end
            end
        end

        % Spike-Contrast implementation for each trialType + cellType pair
        cellTypes = unique(cellTypesVec_spiketrains);
        for cellType_ind = 1:length(cellTypes)
            cellType = cellTypes{cellType_ind};
            respVals = unique(ResponsivityVec_spiketrains);
            for resp_ind = 1:length(respVals)
                StimResp = respVals{resp_ind};
                for trialType_ind = 1:length(trialTagsLabels)
                    trialType = trialTagsLabels{trialType_ind};
    
                    % index the right set of spike trains
                    SpikeTrains = allSpikeTrains(:, strcmp(cellTypesVec_spiketrains,cellType) & strcmp(ResponsivityVec_spiketrains,StimResp) & strcmp(trialTypesVec_spiketrains,trialType));
    
                    [S,PREF] = SpikeContrast(SpikeTrains, recLength, binSizes(1), binSizes(2)); % S is the synchrony index, 0 <= S <= 1
                    s_curve = PREF.s;
                    s_bins = PREF.bins;
            
                    S_vec(end+1,1) = S;
                    s_curves(end+1,:) = s_curve;
                    groupsVec{end+1,1} = groupName;
                    cellTypesVec{end+1,1} = cellType;
                    trialTypesVec{end+1,1} = trialType;
                    responsivityVec{end+1,1} = StimResp;
            
                    % % plot synchrony curve per recording (sanity check)
                    % figure;
                    % plot(s_bins, s_curve);
                    % title(recName);
                end
            end
        end

    end
end

%% Bar plot
figure;

g = gramm('x',trialTypesVec, 'y',S_vec, 'color',groupsVec);
g.facet_grid(cellTypesVec, responsivityVec);
g.stat_summary('type','bootci', 'geom',{'bar','black_errorbar'}, 'setylim',true);
g.set_names('x','', 'y','Synchrony Index', 'Color','', 'Row','', 'Column','');
%g.no_legend;
g.draw();
if plot_points
    g.update('x',trialTypesVec, 'y',S_vec, "color",groupsVec);
    g.geom_jitter();
    g.set_color_options('lightness',40);
    g.set_point_options("markers","^", "base_size",3);
    g.no_legend;
    g.draw;
end

%% Mean synchrony curves
figure;
g = gramm('x',s_bins', 'y',s_curves, 'color',groupsVec);
g.facet_grid(cellTypesVec, responsivityVec);
g.fig(trialTypesVec);
g.stat_summary('type','bootci', 'setylim',true);
g.set_names('x','Time bin (s)', 'y','Synchrony', 'Color','', 'Row','', 'Column','');
% if numel(trialTagsLabels) > 1
%     g.set_order_options('Row',trialTagsLabels);
% end
g.draw;

%% Individual synchrony curves
figure;
g = gramm('x',s_bins, 'y',s_curves, 'color',groupsVec);
g.facet_grid(cellTypesVec, responsivityVec);
g.fig(trialTypesVec);
g.geom_line();
g.set_names('x','Time bin (s)', 'y','Synchrony', 'Color','', 'Row','', 'Column','');
% if numel(trialTagsLabels) > 1
%     g.set_order_options('Row',trialTagsLabels);
% end
g.draw;

%% Save data table
data_table_syn = table(groupsVec, cellTypesVec, responsivityVec, trialTypesVec, S_vec, 'VariableNames',{'Group','CellType','StimResponsivity','TrialType','Synchrony'});

end