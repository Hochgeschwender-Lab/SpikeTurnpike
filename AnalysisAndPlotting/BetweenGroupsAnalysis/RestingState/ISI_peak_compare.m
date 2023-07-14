function stats_out = ISI_peak_compare(all_data, cell_type)
% Plot and compare ISI coefficient of variation (CV) between groups.
%
% INPUTS
% - all_data
% - cell_type: cell type to compare CVs for. e.g., 'MSN'

groupNames = fieldnames(all_data);

groupsVec = {};
ISI_peak_vec = [];

for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};

    mouseNames = fieldnames(all_data.(groupName));

    for mouseNum = 1:length(mouseNames)
        mouseName = mouseNames{mouseNum};

        cellIDs = fieldnames(all_data.(groupName).(mouseName));

        for cellID_num = 1:length(cellIDs)
            cellID = cellIDs{cellID_num};

            thisCellType = all_data.(groupName).(mouseName).(cellID).Cell_Type;
            isSingleUnit = all_data.(groupName).(mouseName).(cellID).IsSingleUnit;
            if strcmp(thisCellType, cell_type) && isSingleUnit
                groupsVec{end+1} = groupName;
                ISI_peak_vec(end+1) = all_data.(groupName).(mouseName).(cellID).ISI_pdf_peak_xy(1);
            end
        end
    end
end

%% Plotting
figure;

g = gramm('x',groupsVec, 'y',ISI_peak_vec, 'color',groupsVec);
g.stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'setylim',true);
g.set_names('x','', 'y','ISI CV', 'Color','');

g.draw();

end

