function FSL_by_group_barplot(all_data,trialTagsLabels,trialType,excludeMUA)
% Bar plot of only positively-modulated units' first spike latencies,
% overlayed between groups.

groupNames = fieldnames(all_data);

trialTypeInd = find(strcmp(trialTagsLabels,trialType));

FSL_vec = [];
FSL_reliability_vec = [];
cellTypesVec = {};
groupVec = {};

for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};
    mouseNames = fieldnames(all_data.(groupName));

    for mouseNum = 1:length(mouseNames)
        mouseName = mouseNames{mouseNum};
        cellIDs = fieldnames(all_data.(groupName).(mouseName));

        for cellID_num = 1:length(cellIDs)
            cellID = cellIDs{cellID_num};

            isSingleUnit = all_data.(groupName).(mouseName).(cellID).IsSingleUnit;
            responsivityNum = all_data.(groupName).(mouseName).(cellID).StimResponsivity;
            if (responsivityNum == 1) && isSingleUnit % unit is positively modulated by stim
                groupVec{end+1,1} = groupName;
                FSL_vec(end+1,1) = all_data.(groupName).(mouseName).(cellID).FirstSpikeLatency(trialTypeInd);
                FSL_reliability_vec(end+1,1) = all_data.(groupName).(mouseName).(cellID).FirstSpikeLatency_Reliability(trialTypeInd);
                cellTypesVec{end+1,1} = all_data.(groupName).(mouseName).(cellID).Cell_Type;
            end
        end
    end
end

%% Plotting
cell_types = unique(cellTypesVec);

figure;

g(1,1) = gramm('x',cellTypesVec, 'y',FSL_vec, 'color',groupVec);
g(1,1).stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'setylim','true');
g(1,1).set_names('x','', 'y','First Spike Latency (ms)');
g(1,1).set_layout_options("legend",false);

g(1,2) = gramm('x',cellTypesVec, 'y',FSL_reliability_vec, 'color',groupVec);
g(1,2).stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'setylim','true');
g(1,2).set_names('x','', 'y','First Spike Latency Reliability', 'Color','');
g(1,2).set_layout_options("legend",false);

g.draw;

end

