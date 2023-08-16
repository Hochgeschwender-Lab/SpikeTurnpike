function StimPhaseLocking(all_data,StimVoltageTrace_ms,minNspikes)
% Computes a spike-triggered average and pairwise phase consistency 2
% (PPC2) for each unit w.r.t. the stimulation voltage trace, not the LFP.


if nargin < 3 || isempty(minNspikes)
    minNspikes = 2;
end

% Extract phase angle of stimulator voltage trace
StimTrace_phase = angle(hilbert(StimVoltageTrace_ms, 2^nextpow2(length(StimVoltageTrace_ms))));
StimTrace_phase = StimTrace_phase(1:length(StimVoltageTrace_ms));
hist_lim = rad2deg(max(abs(StimTrace_phase)));
hist_lims = [-hist_lim hist_lim];

% initialize vectors for plotting/stats, all of length numUnits
PrefPhase_vec = [];
PPC_vec = [];
groupVec = [];
cellTypeVec = {};
responsivityVec = {};
cellTypeVec2 = {};
%trialTypeVec = {};

groupNames = fieldnames(all_data);
for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};

    recNames = fieldnames(all_data.(groupName));

    for recNum = 1:length(recNames)
        recName = recNames{recNum};

        fprintf('  Working on recording %s\n', recName);

        cellIDs = fieldnames(all_data.(groupName).(recName));
        for cellID_num = 1:length(cellIDs)
            cellID = cellIDs{cellID_num};

            thisCellType = all_data.(groupName).(recName).(cellID).Cell_Type;
            ISI_violations_percent = all_data.(groupName).(recName).(cellID).ISI_violations_percent;

            if ISI_violations_percent < 1.5
                % get spike trains (binarized) with millisecond indices
                SpikeTrains_ms = all_data.(groupName).(recName).(cellID).SpikeTrains_stim_ms'; % [time x trials] binarized spike data

                % only proceed if the unit spikes at least two
                % times in the stimulation period
                if sum(SpikeTrains_ms,'all') >= minNspikes
                    StimTrace_phase_i = repmat(StimTrace_phase, 1, size(SpikeTrains_ms,2));
                    StimTrace_phase_i(~SpikeTrains_ms) = 0;
                    [~, spikeTrialNums, spikePhases] = find(StimTrace_phase_i);

                    % calculate preferred phase as the circular
                    % mean of all spike phases
                    PrefPhase_vec(end+1,1) = circ_mean(spikePhases);

                    % calculate Pairwise Phase Consistency P^_2
                    % as described in Vinck et al. (2012, J Comp
                    % Neurosci)
                    PPC_vec(end+1,1) = PPC2(spikePhases, spikeTrialNums);
                    clear spikePhases spikeTrialNums

                    % populate other vectors for plotting
                    groupVec{end+1,1} = groupName;
                    cell_type = all_data.(groupName).(recName).(cellID).Cell_Type;
                    cellTypeVec{end+1,1} = cell_type;
                    responsivityNum = all_data.(groupName).(recName).(cellID).StimResponsivity;
                    if responsivityNum == 1
                        responsivityVec{end+1,1} = '+';
                        cellTypeVec2{end+1,1} = strcat('+',cell_type);
                    elseif responsivityNum == 0
                        responsivityVec{end+1,1} = 'nr';
                        cellTypeVec2{end+1,1} = strcat('ns',cell_type);
                    else
                        responsivityVec{end+1,1} = '-';
                        cellTypeVec2{end+1,1} = strcat('-',cell_type);
                    end

                end
            end
        end
    end
end

%% PPC2 bar plot
figure;
g = gramm('x',groupVec, 'y',PPC_vec, 'color',groupVec);
g.facet_grid(cellTypeVec, responsivityVec, "scale","independent");
g.stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'setylim',true);
g.set_names('x','', 'y','Spike-Stim PPC2', 'color','', 'column','')
g.no_legend;
g.draw;

g.update('x',groupVec, 'y',PPC_vec, 'color',groupVec);
g.geom_jitter();
g.set_color_options('lightness',40);
g.set_point_options('base_size',2.5);
g.no_legend;
g.draw;


%% Plot preferred phase histogram
figure;
g = gramm('x',PrefPhase_vec, 'color',groupVec);
g.facet_grid(cellTypeVec, responsivityVec, 'scale','independent');
% g.stat_bin('geom','overlaid_bar', 'normalization','pdf', 'edges',-180:20:180);
g.stat_bin('geom','overlaid_bar', 'normalization','count', 'edges',-hist_lim:1:hist_lim);
%g.set_polar("closed",0);
g.set_order_options("column",{'+','nr','-'});
g.set_names('x','Preferred Phase', 'y','PDF', 'color','', 'row','', 'column','');
g.draw;

end