function plot_whole_rec_rasters(all_data, cell_types, sortByVar, recs)
% For each recording, generates a raster plot with units along the y axis
% color-coded by cell type.
%
% INPUTS
% - all_data
% - cell_types: cell array of cell types to use for plots
% - sortByVar: string, variable in all_data to sort units by along the y
%              axis. If empty, will not sort by a variable.
% - recs: names of recordings to plot. If empty, only plots all recordings.

groupNames = fieldnames(all_data);

% arrays to populate
spike_trains = {};
cellTypes_vec = {};
recNames_vec = {};
groups_vec = {};
sortByVar_vec = [];
cellID_vec = {};

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

                cellTypes_vec{end+1,1} = thisCellType;
                spike_trains{end+1,1} = spikeTimes_s;
                recNames_vec{end+1,1} = recName;
                groups_vec{end+1,1} = groupName;
                cellID_vec{end+1,1} = cellID;

                if ~isempty(sortByVar)
                    sortByVar_vec(end+1) = all_data.(groupName).(recName).(cellID).(sortByVar);
                    % BurstsInfo = all_data.(groupName).(recName).(cellID).BurstsInfo;
                    % sortByVar_vec(end+1) = BurstsInfo.proportion_spikes_in_bursts;
                end
            end
        end
    end
end

%% Plot big rasters
figure;
g = gramm('x',spike_trains, 'color',cellTypes_vec, 'lightness',sortByVar_vec, 'label',cellID_vec);
%g = gramm('x',spike_trains, 'lightness',sortByVar_vec, 'label',cellID_vec);
g.facet_grid(recNames_vec,[], 'scale','free');
g.fig(groups_vec);
g.geom_raster();
%g.set_order_options('color',cell_types);
g.set_names('x','Time (s)', 'y','Unit', 'color','', 'row','', 'fig','Group');

g.set_text_options('interpreter','none');

% if ~isempty(sortByVar)
%     g.set_order_options("x",sortByVar_vec);
% end

g.draw;

%% Plot with only specified recordings
recInds_matrix = [];
for rec_ind = 1:length(recs)
    rec = recs{rec_ind};
    recInds_matrix(end+1,:) = strcmp(recNames_vec, rec);
end
recInds = find(sum(recInds_matrix,1));

figure;
%g = gramm('x',spike_trains(recInds), 'color',cellTypes_vec(recInds), 'lightness',sortByVar_vec(recInds));
%g = gramm('x',spike_trains(recInds), 'color',groups_vec(recInds), 'lightness',sortByVar_vec(recInds));
g = gramm('x',spike_trains(recInds), 'color',groups_vec(recInds));

g.facet_grid(recNames_vec(recInds),[], 'scale','free');
g.geom_raster();
%g.set_order_options('color',cell_types);
g.set_names('x','Time (s)', 'y','Unit', 'color','', 'row','');

g.set_text_options('interpreter','none');

g.axe_property('XLim', [1000 1200]);

g.draw;

end

