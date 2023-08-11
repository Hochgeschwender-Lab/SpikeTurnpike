function PlotDCW(MonoConnectionsTable)
% Plot directed connection weight (DCW), see Jia et al. (2021).

groups = unique(MonoConnectionsTable.group);
cellTypePairs = unique(MonoConnectionsTable.CellTypes);

% initialize vectors for plotting
dcw_vec = [];
groups_vec = {};
cellTypePairs_vec = {};
EorI_vec = {};

for groupNum = 1:length(groups)
    for cellTypePairNum = 1:length(cellTypePairs)
        inds = ismember(MonoConnectionsTable.group,groups{groupNum}) & ismember(MonoConnectionsTable.CellTypes,cellTypePairs{cellTypePairNum});

        dcw_vec = [dcw_vec; MonoConnectionsTable.DCW(inds)];
        %[groups_vec{end+1:end+1+sum(inds)}] = MonoConnectionsTable.group{inds};
        groups_vec = [groups_vec; MonoConnectionsTable.group(inds)];
        cellTypePairs_vec = [cellTypePairs_vec; MonoConnectionsTable.CellTypes(inds)];
        EorI_vec = [EorI_vec; MonoConnectionsTable.EorI(inds)];
    end
end

%% Plotting
figure;
g = gramm('x',EorI_vec, 'y',dcw_vec, 'color',groups_vec, 'column',cellTypePairs_vec);
g.stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'setylim',true);
g.set_names('x','', 'y','Directed Connection Weight', 'color','', 'column','')
g.no_legend;
g.draw;

g.update('x',EorI_vec, 'y',dcw_vec, 'color',groups_vec, 'column',cellTypePairs_vec);
g.geom_point();
g.set_color_options('lightness',40);
g.set_point_options('base_size',2.5);
% g.no_legend;
g.draw;

end

