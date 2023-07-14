function electrode_spike_identification = spike_identification_fun(filtered_electrode_data, file_directory)

% Get the number of rows and columns of the electrode_data
[row_filtered, column_filtered] = size(filtered_electrode_data);

%% Identification of spikes
standard_deviation_array = zeros(row_filtered, 1);

for n_electrodes=1:row_filtered
    median_electrode_data = median(filtered_electrode_data(n_electrodes, :));
    absolute_electrode_data_minus_median = abs(filtered_electrode_data(n_electrodes, :) - median_electrode_data);
    
    standard_deviation_array(n_electrodes) = median(absolute_electrode_data_minus_median/0.6745);
end

electrode_spike_identification = zeros(row_filtered, column_filtered);
for n_electrodes=1:row_filtered
    spike_criteria = - 3 * standard_deviation_array(n_electrodes);
    spike_or_not_spike_array = filtered_electrode_data(n_electrodes, :) < spike_criteria;
    
    for n_samples=1:column_filtered
        electrode_spike_identification(n_electrodes, n_samples) = double(spike_or_not_spike_array(1, n_samples));
    end
    fprintf('Finished spike identification Electrode %d\n', n_electrodes)
end

 save(strcat(file_directory, 'identified_electrode_spike.mat'),...
     'electrode_spike_identification', '-v7.3')

end

