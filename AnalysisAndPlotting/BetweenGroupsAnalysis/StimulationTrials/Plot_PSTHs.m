function Plot_PSTHs(all_data,relative_time_ms,trialTagsLabels,PSTH_type,excludeMUA)
%
% PSTH_type: 'raw' or 'conv'
% excludeMUA: 0 or 1

groupNames = fieldnames(all_data);

% all of these should have n_rows = n_units
all_PSTHs = []; % [units x time] matrix of trial-average PSTHs
groupsVec = {};
trialTagsVec = {};
cellTypesVec = {}; % elements are e.g. 'FS', 'RS', ...
cellTypesVec2 = {}; % elements are e.g. '+FS', '-FS', '+RS', 'nsRS', ...
responsivityVec = {};
layerVec = {};

for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};
    mouseNames = fieldnames(all_data.(groupName));

    for mouseNum = 1:length(mouseNames)
        mouseName = mouseNames{mouseNum};
        cellIDs = fieldnames(all_data.(groupName).(mouseName));

        for cellID_num = 1:length(cellIDs)
            cellID = cellIDs{cellID_num};

            IsSingleUnit = all_data.(groupName).(mouseName).(cellID).IsSingleUnit;
            ISI_violations_percent = all_data.(groupName).(mouseName).(cellID).ISI_violations_percent;
            
            %if (excludeMUA && IsSingleUnit) || ~excludeMUA
            if (excludeMUA && (ISI_violations_percent <= 1)) || ~excludeMUA
                for trialTagInd = 1:length(trialTagsLabels)
                    if strcmp(PSTH_type, 'raw')
                        all_PSTHs(end+1,:) = all_data.(groupName).(mouseName).(cellID).PSTHs_raw(trialTagInd,:);
                    else
                        all_PSTHs(end+1,:) = all_data.(groupName).(mouseName).(cellID).PSTHs_conv(trialTagInd,:);
                    end
                    groupsVec{end+1,1} = groupName;
    
                    cell_type = all_data.(groupName).(mouseName).(cellID).Cell_Type;
                    cellTypesVec{end+1,1} = cell_type;
                    trialTagsVec{end+1,1} = trialTagsLabels{trialTagInd};
                    layerVec{end+1,1} = all_data.(groupName).(mouseName).(cellID).LaminarLabel;
        
                    responsivityNum = all_data.(groupName).(mouseName).(cellID).StimResponsivity;
                    if responsivityNum == 1
                        responsivityLabel = '+';
                        responsivityVec{end+1,1} = responsivityLabel;
                    elseif responsivityNum == 0
                        responsivityLabel = 'nr';
                        responsivityVec{end+1,1} = responsivityLabel;
                    else
                        responsivityLabel = '-';
                        responsivityVec{end+1,1} = responsivityLabel;
                    end
    
                    cellTypesVec2{end+1,1} = strcat(responsivityLabel,cell_type);
                end
            end

        end
    end
end

%% Plotting
figure;
g = gramm('x',relative_time_ms, 'y',all_PSTHs, 'color',groupsVec);
%g.facet_grid(cellTypesVec, responsivityVec, "scale","free_y");
g.facet_grid(layerVec, cellTypesVec2, "scale","free_y");
g.fig(trialTagsVec);
g.stat_summary('type','sem', 'setylim','true');
%g.geom_line();
g.set_names('x','Time (ms)', 'y','Instantaneous FR (Hz)', 'Color','', 'Row','', 'Column','');
g.set_order_options('row',{'SG','L4','IG'});
g.draw;

%% TO DO: Print sample sizes as a table (cell type x modulation)

end

