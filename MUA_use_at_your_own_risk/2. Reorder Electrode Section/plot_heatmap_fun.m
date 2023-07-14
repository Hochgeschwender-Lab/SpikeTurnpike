function plot_heatmap_fun(density_matrix, particular_file_name, file_directory_to_save_data)

%{
    Plot the density_matrix data binned to 1 second intervals as a heat map
    to observe the entire electrode data as a whole. These figures will be
    stored inside the sub-folder "AllElectrodes" within the folder 
    "figures".
%}

% Save all the figures to a separate folder called "figures/AllElectrodes".
check_if_folder_exists(strcat(file_directory_to_save_data, 'figures/AllElectrodes/'));

% Get the dimension information of the spike_density variable
[row_density_matrix, column_density_matrix] = size(density_matrix);

% Total time-span of data
start_time = 0;
time_span = column_density_matrix; 

% Time bins every how many seconds?
time_bins = 1; % bins spike values every "time_bins" second

time = start_time:time_bins:time_span; % s

electrodes = linspace(1, row_density_matrix, row_density_matrix);
figure('visible', 'off');
imagesc(time, electrodes, density_matrix);
colorbar;
filename = strcat(file_directory_to_save_data, 'figures/AllElectrodes/', particular_file_name, '.png');
saveas(gcf, filename);

end

