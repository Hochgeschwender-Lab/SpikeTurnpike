function PlotMeanCCGs(MonoConnectionsTable,MonoConnCCGs,t)
%PLOTMONOSYNAPTICCONNECTIONS Summary of this function goes here
%   Detailed explanation goes here

% groups = unique(MonoConnectionsTable.group);
% cellTypePairs = unique(MonoConnectionsTable.CellTypes);
% 
% for cellTypePairNum = 1:length(cellTypePairs)
%     figure;
%     T = tiledlayout(1,2, "TileSpacing","compact");
%     title(T, cellTypePairs{cellTypePairNum});
%     maxmaxZ = 0;
%     minminZ = 0;
%     for groupNum = 1:length(groups)
%         nexttile(T);
% 
%         inds_E = ismember(MonoConnectionsTable.group,groups{groupNum}) & ismember(MonoConnectionsTable.CellTypes,cellTypePairs{cellTypePairNum}) & ismember(MonoConnectionsTable.EorI,'E');
%         CCGsToPlot_E = MonoConnCCGs(:,inds_E);
% 
%         inds_I = ismember(MonoConnectionsTable.group,groups{groupNum}) & ismember(MonoConnectionsTable.CellTypes,cellTypePairs{cellTypePairNum}) & ismember(MonoConnectionsTable.EorI,'I');
%         CCGsToPlot_I = MonoConnCCGs(:,inds_I);
% 
%         [~, maxInds] = max(CCGsToPlot_E(t>0 & t<=10, :),[],1);
%         [~, I_sort] = sort(maxInds);
%         CCGsToPlot_E = CCGsToPlot_E(:,I_sort);
% 
%         [~, minInds] = min(CCGsToPlot_I(t>0 & t<=10, :),[],1);
%         [~, I_sort] = sort(minInds);
%         CCGsToPlot_I = CCGsToPlot_I(:,I_sort);
% 
%         CCGsToPlot = zscore([CCGsToPlot_E CCGsToPlot_I]);
% 
%     end
% end

peaks_vec = max(abs(MonoConnCCGs(t>0 & t<=5, :)),[],1);

%% Plotting
% overlaid mean CCGs
figure;
g = gramm('x',t, 'y',MonoConnCCGs', 'color',MonoConnectionsTable.group);
g.facet_grid(MonoConnectionsTable.EorI, MonoConnectionsTable.CellTypes, 'scale','independent');
g.stat_summary('type','sem', 'setylim',true);
%g.geom_line;
g.draw;


% all CCGs
figure;
g = gramm('x',t, 'y',MonoConnCCGs', 'color',MonoConnectionsTable.group);
g.facet_grid(MonoConnectionsTable.EorI, MonoConnectionsTable.CellTypes, 'scale','independent');
%g.stat_summary('type','sem', 'setylim',true);
g.geom_line;
g.draw;


% bar plot of peak values
figure;
g = gramm('x',MonoConnectionsTable.group, 'y',peaks_vec, 'color',MonoConnectionsTable.group);
g.facet_grid(MonoConnectionsTable.EorI, MonoConnectionsTable.CellTypes, 'scale','independent');
g.stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'setylim',true);
g.set_names('x','', 'y','CCG Peak (1 ms <= t <= 5 ms)', 'color','', 'column','')
g.set_title('Monosynaptic CCG Peak');
g.no_legend;
g.draw;

g.update('x',MonoConnectionsTable.group, 'y',peaks_vec, 'color',MonoConnectionsTable.group);
g.geom_jitter();
g.set_color_options('lightness',40);
g.set_point_options('base_size',2.5);
g.no_legend;
g.draw;


% bar plot of excess synchrony for significantly synchronous pairs
figure;
g = gramm('x',MonoConnectionsTable.group, 'y',MonoConnectionsTable.ExcessSynchrony, 'color',MonoConnectionsTable.group, 'subset',strcmp(MonoConnectionsTable.EorI,'S'));
g.facet_grid([], MonoConnectionsTable.CellTypes, 'scale','independent');
%g.stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'setylim',true);
g.geom_jitter();
g.set_names('x','', 'y','Excess Synchrony', 'color','', 'column','')
g.no_legend;
g.set_title('Excess Synchrony: Significant Synchronous Pairs');
g.draw;

% g = gramm('x',MonoConnectionsTable.group, 'y',MonoConnectionsTable.ExcessSynchrony, 'color',MonoConnectionsTable.group, 'subset',strcmp(MonoConnectionsTable.EorI,'S'));
% g.geom_jitter();
% g.set_color_options('lightness',40);
% g.set_point_options('base_size',2.5);
% g.no_legend;
% g.draw;


% bar plot of excess synchrony for all pairs
figure;
g = gramm('x',MonoConnectionsTable.layers, 'y',MonoConnectionsTable.ExcessSynchrony, 'color',MonoConnectionsTable.group);
g.facet_grid([], MonoConnectionsTable.CellTypes, 'scale','independent');
g.stat_summary('type','bootci', 'geom',{'bar','black_errorbar'}, 'setylim',true);
% g.geom_jitter();
g.set_names('x','', 'y','Excess Synchrony', 'color','', 'column','')
g.no_legend;
g.set_title('Excess Synchrony: All Pairs');
g.draw;

% g = gramm('x',MonoConnectionsTable.group, 'y',MonoConnectionsTable.ExcessSynchrony, 'color',MonoConnectionsTable.group);
% g.geom_jitter();
% g.set_color_options('lightness',40);
% g.set_point_options('base_size',2.5);
% g.no_legend;
% g.draw;

end