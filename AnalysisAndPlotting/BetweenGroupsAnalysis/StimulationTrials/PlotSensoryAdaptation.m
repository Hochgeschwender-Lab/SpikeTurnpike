function PlotSensoryAdaptation(all_data,trialTagsLabels,userGroupNames,groupsPlottingOrder,trialTagToPlot)
% Plots the response ratio for each deflection. Only uses units that are
% positively modulated by the stimulus and pass the ISI violations quality
% control check (no MUA).
%
% INPUTS
%   userGroupNames: {'group1','group2',...} to rename groups in the same
%       order that they appear in the all_data struct.
%   groupsPlottingOrder: {'group1','group2',...}, with the names given in
%       userGroupNames, to determine the order of plotting in Gramm.

groupNames = fieldnames(all_data);

trialTagInd = find(strcmp(trialTagsLabels, trialTagToPlot));

groupsVec = {};
cellTypesVec = {};
responseRatiosVec = [];
deflectionNumVec = [];

for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};
    userGroupName = userGroupNames{groupNum};
    mouseNames = fieldnames(all_data.(groupName));

    for mouseNum = 1:length(mouseNames)
        mouseName = mouseNames{mouseNum};
        cellIDs = fieldnames(all_data.(groupName).(mouseName));

        for cellID_num = 1:length(cellIDs)
            cellID = cellIDs{cellID_num};

            StimResponsivity = all_data.(groupName).(mouseName).(cellID).StimResponsivity;
            IsSingleUnit = all_data.(groupName).(mouseName).(cellID).IsSingleUnit;
            if IsSingleUnit && (StimResponsivity == 1)
                responseRatios = all_data.(groupName).(mouseName).(cellID).DeflectionResponseRatios(trialTagInd,:);
                for deflectionNum = 1:length(responseRatios)
                    responseRatiosVec(end+1,1) = responseRatios(deflectionNum);
                    deflectionNumVec(end+1,1) = deflectionNum;
                    groupsVec{end+1,1} = userGroupName;
                    cellTypesVec{end+1,1} = all_data.(groupName).(mouseName).(cellID).Cell_Type;
                end
            end
        end
    end
end

figure;

g = gramm("x",deflectionNumVec, "y",responseRatiosVec, "color",groupsVec, "column",cellTypesVec);
g.stat_summary('type','sem', 'geom',{'line','black_errorbar','point'}, 'setylim',true);
%g.geom_point();
g.set_names('x','Deflection', 'y','Response Ratio', 'Color','', 'Column','', 'Row','');
g.set_order_options("color",groupsPlottingOrder);

g.draw();

end