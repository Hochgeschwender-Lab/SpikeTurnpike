function PlotMonosynapticConnections(MonoConnectionsTable,MonoConnCCGs,t)
%PLOTMONOSYNAPTICCONNECTIONS Summary of this function goes here
%   Detailed explanation goes here

groups = unique(MonoConnectionsTable.group);
cellTypePairs = unique(MonoConnectionsTable.CellTypes);

for cellTypePairNum = 1:length(cellTypePairs)
    figure;
    T = tiledlayout(1,2, "TileSpacing","compact");
    title(T, cellTypePairs{cellTypePairNum});
    maxmaxZ = 0;
    minminZ = 0;
    for groupNum = 1:length(groups)
        nexttile(T);
        
        inds_E = ismember(MonoConnectionsTable.group,groups{groupNum}) & ismember(MonoConnectionsTable.CellTypes,cellTypePairs{cellTypePairNum}) & ismember(MonoConnectionsTable.EorI,'E');
        CCGsToPlot_E = MonoConnCCGs(:,inds_E);

        inds_I = ismember(MonoConnectionsTable.group,groups{groupNum}) & ismember(MonoConnectionsTable.CellTypes,cellTypePairs{cellTypePairNum}) & ismember(MonoConnectionsTable.EorI,'I');
        CCGsToPlot_I = MonoConnCCGs(:,inds_I);
        
        [~, maxInds] = max(CCGsToPlot_E(t>0 & t<=10, :),[],1);
        [~, I_sort] = sort(maxInds);
        CCGsToPlot_E = CCGsToPlot_E(:,I_sort);

        [~, minInds] = min(CCGsToPlot_I(t>0 & t<=10, :),[],1);
        [~, I_sort] = sort(minInds);
        CCGsToPlot_I = CCGsToPlot_I(:,I_sort);
        
        CCGsToPlot = zscore([CCGsToPlot_E CCGsToPlot_I]);

        maxZ = max(CCGsToPlot,[],'all');
        if maxZ > maxmaxZ
            maxmaxZ = maxZ;
        end
        minZ = min(CCGsToPlot,[],'all');
        if minZ < minminZ
            minminZ = minZ;
        end

        plot_matrix(CCGsToPlot, t, 1:size(CCGsToPlot,2), 'n');
        set(gca, 'YDir','reverse');
        title(gca, groups{groupNum});
        colorbar(gca,'off');
        xlabel('');
        ylabel('');
        set(gca,'YTickLabel',[]);

        xline([-10 10],'r');
        xline(0);
    end

    %% create shared clims and colorbar
    for groupNum = 1:length(groups)
        nexttile(T,groupNum);
        clim(gca,[minminZ maxmaxZ]);
    end
    c = colorbar(gca, "eastoutside", "Limits",[minminZ maxmaxZ]);
    c.Label.String = "Z-Score";
    xlabel(T,'Lag (ms)');
    ylabel(T,'Monosynaptically Connected Pairs');
end

end