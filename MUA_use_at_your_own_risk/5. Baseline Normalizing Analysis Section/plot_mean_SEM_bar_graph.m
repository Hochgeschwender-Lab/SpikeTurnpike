function plot_mean_SEM_bar_graph(mean_extracted_data, SEM_extracted_data, ...
    electrodes_order, n_electrodes, total_time_range, key_word, file_directory)

new_electrode_number = n_electrodes;
original_electrode_number = electrodes_order(new_electrode_number);

% Save all the figures to a separate folder called "figures/SEMBarGraph".
check_if_folder_exists(strcat(file_directory, 'figures/SEM_BarGraph/', key_word));

time_stamp_bar_graph = 1:total_time_range;
NOT_normalized_error_high = SEM_extracted_data;
NOT_normalized_error_low = -SEM_extracted_data;

bar(time_stamp_bar_graph, mean_extracted_data);
hold on;
errorbar(time_stamp_bar_graph, mean_extracted_data, NOT_normalized_error_low, NOT_normalized_error_high);
xlabel('Electrode ' + string(original_electrode_number))
hold off;
filename = strcat(file_directory, 'figures/SEM_BarGraph/', key_word, '/Electrode_', string(original_electrode_number), '_SEMBarGraph.png');
saveas(gcf, filename);

end

