function Plot_PSTHs(all_data,relative_time_ms,trialTagsLabels,trialTagToPlot,PSTH_type,excludeMUA)
%
% PSTH_type: 'raw' or 'conv'
% excludeMUA: 0 or 1

groupNames = fieldnames(all_data);

% all of these should have n_rows = n_units
all_PSTHs = []; % [units x time] matrix of trial-average PSTHs
groupsVec = {};
cellTypesVec = {};
responsivityVec = {};

% To get the PSTH for the right type of trial
trialTagInd = find(strcmp(trialTagsLabels, trialTagToPlot));

if isempty(trialTagInd)
    error("Invalid trialTagToPlot. Check your trialTagsLabels for the available options and try again.")
end

for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};
    mouseNames = fieldnames(all_data.(groupName));

    for mouseNum = 1:length(mouseNames)
        mouseName = mouseNames{mouseNum};
        cellIDs = fieldnames(all_data.(groupName).(mouseName));

        for cellID_num = 1:length(cellIDs)
            cellID = cellIDs{cellID_num};

            IsSingleUnit = all_data.(groupName).(mouseName).(cellID).IsSingleUnit;
            
            if (excludeMUA && IsSingleUnit) || ~excludeMUA
                if strcmp(PSTH_type, 'raw')
                    all_PSTHs(end+1,:) = all_data.(groupName).(mouseName).(cellID).PSTHs_raw(trialTagInd,:);
                else
                    all_PSTHs(end+1,:) = all_data.(groupName).(mouseName).(cellID).PSTHs_conv(trialTagInd,:);
                end
                groupsVec{end+1,1} = groupName;
                cellTypesVec{end+1,1} = all_data.(groupName).(mouseName).(cellID).Cell_Type;
    
                responsivityNum = all_data.(groupName).(mouseName).(cellID).StimResponsivity;
                if responsivityNum == 1
                    responsivityVec{end+1,1} = '+';
                elseif responsivityNum == 0
                    responsivityVec{end+1,1} = 'NR';
                else
                    responsivityVec{end+1,1} = '-';
                end
            end

        end
    end
end

%% Plotting
figure;
g = gramm('x',relative_time_ms, 'y',all_PSTHs, 'color',groupsVec);
g.facet_grid(cellTypesVec, responsivityVec, "scale","free_y");
g.stat_summary('type','sem', 'setylim','true');
%g.geom_line();
g.set_names('x','Time (ms)', 'y','Instantaneous FR (Hz)', 'Color','', 'Row','', 'Column','');
%g.set_order_options("row",['+','NR','-']);
g.draw;

%% TO DO: Print sample sizes as a table (cell type x modulation)

end

