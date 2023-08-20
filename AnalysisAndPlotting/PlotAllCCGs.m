function PlotAllCCGs(MonoConnectionsTable, MonoConnCCGs, SimultaneousABs, PointwiseABs, t)
%PLOTALLSIGNIFICANTCCGS Summary of this function goes here
%   Detailed explanation goes here

groupNames = unique(MonoConnectionsTable.group);
for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};
    dataTable_group = MonoConnectionsTable(strcmp(MonoConnectionsTable.group, groupName),:);
    CCGs_group = MonoConnCCGs(:,strcmp(MonoConnectionsTable.group, groupName));
    sim_ABs_group = SimultaneousABs(:,strcmp(MonoConnectionsTable.group, groupName),:);
    pt_ABs_group = PointwiseABs(:,strcmp(MonoConnectionsTable.group, groupName),:);
    
    recNames = unique(dataTable_group.recording);
    for recNum = 1:length(recNames)
        recName = recNames{recNum};
        dataTable_rec = dataTable_group(strcmp(dataTable_group.recording, recName),:);
        CCGs_rec = CCGs_group(:,strcmp(dataTable_group.recording, recName),:);
        sim_ABs_rec = sim_ABs_group(:,strcmp(dataTable_group.recording, recName),:);
        pt_ABs_rec = pt_ABs_group(:,strcmp(dataTable_group.recording, recName),:);

        figure('Name',recName);
        tiledlayout('flow');
        
        numPairs = height(dataTable_rec);
        for pair_ind = 1:numPairs
            if dataTable_rec.Significance(pair_ind) == 1
                % plot the CCG
                ccg_i = squeeze(CCGs_rec(:,pair_ind,:));
                AB_sim = squeeze(sim_ABs_rec(:,pair_ind,:));
                AB_pt = squeeze(pt_ABs_rec(:,pair_ind,:));

                cellID_A = dataTable_rec.UnitA(pair_ind);
                cellID_B = dataTable_rec.UnitB(pair_ind);

                cellTypes = dataTable_rec.CellTypes(pair_ind);

                EorI = dataTable_rec.EorI(pair_ind);

                nexttile;
                plot(t,AB_sim,"Color","#D95319"); hold on; plot(t,AB_pt,"Color","#EDB120"); plot(t,ccg_i,"Color","#0072BD"); xline([-5 5]); hold off;
    %                            title(gca,strcat(cellID_A,' -> ',cellID_B,' | ',cell_type_A,' -> ',cell_type_B,' | ',layer_A,'->',layer_B));
                title(gca,strcat(EorI,' | ',cellID_A,'->',cellID_B,' | ',cellTypes));
            end
        end
    end
end

