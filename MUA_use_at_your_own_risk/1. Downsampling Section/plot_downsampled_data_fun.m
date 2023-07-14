function plot_downsampled_data_fun(channel_data_original, channel_data_downsampled, n_electrodes, file_directory)

%{
    Plot the figures separately to observe how the signals are downsampled
    without losing information. These figures will be stored inside the
    sub-folder "DownsampledData" within the folder "figures".

%}

% Save all the figures to a separate folder called "figures/DownsampledData".
check_if_folder_exists(strcat(file_directory, 'figures/DownsampledData/Electrode', string(n_electrodes)));

% Plot the original data before downsampling.
figure('visible', 'off');
plot(channel_data_original(n_electrodes, :))
filename = strcat(file_directory, 'figures/DownsampledData/Electrode', string(n_electrodes), '/Original_Electrode_Data.png');
xlabel('Electrode ' + string(n_electrodes) + ' Original')
ylabel('Recordings')
saveas(gcf, filename); 

% Plot the downsampled data.
figure('visible', 'off');
plot(channel_data_downsampled(n_electrodes, :))
filename = strcat(file_directory, 'figures/DownsampledData/Electrode', string(n_electrodes), '/Downsampled_Electrode_Data.png');
xlabel('Electrode ' + string(n_electrodes) + ' Downsampled')
ylabel('Recordings')
saveas(gcf, filename);

end

