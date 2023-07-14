function PlotAmplitudeByDepth(all_data, channel_positions, electrodes_order)
%PLOTUNITSBYDEPTH Summary of this function goes here
%   Detailed explanation goes here

groupNames = fieldnames(all_data);
yticks_values = sort(channel_positions(:,2));
y_limits = [min(yticks_values), 25];
yticks_labels = flip(electrodes_order);

for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};
    mouseNames = fieldnames(all_data.(groupName));
    for mouseNum = 1:length(mouseNames)
        mouseName = mouseNames{mouseNum};
        cellIDs = fieldnames(all_data.(groupName).(mouseName));

        amplitudes = zeros(length(cellIDs),1);
        depths = zeros(length(cellIDs),1);
        cellTypes = {};

        for cellID_num = 1:length(cellIDs)
            cellID = cellIDs{cellID_num};
            amplitudes(cellID_num) = all_data.(groupName).(mouseName).(cellID).Amplitude;
            depths(cellID_num) = all_data.(groupName).(mouseName).(cellID).Template_Channel_Position(2);
            cellTypes{cellID_num,1} = all_data.(groupName).(mouseName).(cellID).Cell_Type;
        end

        fig = figure();
        fig.Position(3:4) = [400 600];
        gscatter(amplitudes,depths,cellTypes);
        title(gca,strcat(groupName,", ",mouseName), 'Interpreter','none');
        xlabel("Amplitude (uV)");
        ylabel("Channel ID by Depth");
        ylim(y_limits);
        yticks(yticks_values);
        yticklabels(yticks_labels);
        set(gca, 'YGrid', 'on', 'XGrid', 'off');
    end
end

