function plot_filtered_data_fun(filtered_electrode_signal, electrodes_ordered_downsampled_data, original_electrode_number, file_directory_to_save_data)

%{
    Plot the filtered data figures separately to observe how the signals
    are filtered. These figures will be stored inside the
    sub-folder "FilteredData" within the folder "figures".
%}

% Save all the figures to a separate folder called "figures/DownsampledData".
check_if_folder_exists(strcat(file_directory_to_save_data, 'figures/FilteredData/Electrode ', string(original_electrode_number)));

% Plot the data before the application of filter
figure('visible', 'off');
plot(electrodes_ordered_downsampled_data)
xlabel('Original Electrode ' + string(original_electrode_number) + ' NOT filtered')
ylabel('Recordings')
filename = strcat(file_directory_to_save_data, 'figures/FilteredData/Electrode', string(original_electrode_number), '/Before_Filter.png');
saveas(gcf, filename);

% Plot the data after the application of filter
figure('visible', 'off');
plot(filtered_electrode_signal)
xlabel('Original Electrode ' + string(original_electrode_number) + ' FILTERED')
ylabel('Recordings')
filename = strcat(file_directory_to_save_data, 'figures/FilteredData/Electrode', string(original_electrode_number), '/After_Filter.png');
saveas(gcf, filename);

%{
    %{
        This part of the code basically plots a section of the before and
        after filter applied data for visualization and analysis.
    %}

    % Plot a fraction of the data
    lower_range = 2000000;
    upper_range = lower_range + 15000;
    figure;
    plot(electrodes_ordered_downsampled_data(sample_number, lower_range:upper_range))

    figure;
    plot(filtered_electrode_signal(sample_number, lower_range:upper_range))
%}

end

