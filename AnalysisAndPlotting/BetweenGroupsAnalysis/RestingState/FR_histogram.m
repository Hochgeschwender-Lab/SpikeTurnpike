function FRs_table = FR_histogram(all_data, cell_types)
%FR_HISTOGRAM Summary of this function goes here
%   Detailed explanation goes here

groupNames = fieldnames(all_data);

groupsVec = {};
recordingsVec = {};
cellTypesVec = {};
cellIDs_vec = {};
FRs_vec = [];

for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};

    recNames = fieldnames(all_data.(groupName));

    for recNum = 1:length(recNames)
        recName = recNames{recNum};

        cellIDs = fieldnames(all_data.(groupName).(recName));

        for cellID_num = 1:length(cellIDs)
            cellID = cellIDs{cellID_num};

            thisCellType = all_data.(groupName).(recName).(cellID).Cell_Type;
            ISI_violations_percent = all_data.(groupName).(recName).(cellID).ISI_violations_percent;
            
            if any(strcmp(cell_types, thisCellType)) && (ISI_violations_percent <= 1.5)
                groupsVec{end+1,1} = groupName;
                recordingsVec{end+1,1} = recName;
                cellTypesVec{end+1,1} = thisCellType;
                cellIDs_vec{end+1,1} = cellID;

                this_FR = length(all_data.(groupName).(recName).(cellID).SpikeTimes_all) / all_data.(groupName).(recName).(cellID).Recording_Duration;

                FRs_vec(end+1,1) = this_FR;
            end
        end
    end
end

%% plot
figure;
g = gramm('x',FRs_vec, 'color',groupsVec);
g.facet_grid(cellTypesVec, []);
% g.stat_bin('geom','overlaid_bar', 'normalization','pdf', 'edges',-180:20:180);
%g.stat_density("npoints",1000);
g.stat_bin("geom","overlaid_bar", "normalization","count");
%g.set_polar("closed",0);
%g.set_order_options("column",{'+','nr','-'});
g.set_names('x','Firing Rate (Hz)', 'y','Number of Units', 'color','', 'row','', 'column','');
g.draw;

%% make table
varNames = {'group','recording','cell_type','cellID','FR'};
FRs_table = table(groupsVec,recordingsVec,cellTypesVec,cellIDs_vec,FRs_vec, 'VariableNames',varNames);

end