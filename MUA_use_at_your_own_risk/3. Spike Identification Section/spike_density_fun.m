function spike_density = spike_density_fun(electrode_spike_identification, sampling_frequency, bin_spikes_every_x_millisecond, file_directory_to_save_data)

[row_spike_identification, column_spike_identification] = size(electrode_spike_identification);

% Total time-span of data
start_time = 1;
total_experiment_time_s = column_spike_identification/sampling_frequency; % s
total_experiment_time_ms = total_experiment_time_s * 1e3;

after_binning_data_size = length(start_time:bin_spikes_every_x_millisecond:total_experiment_time_ms);

% Conversion from milliseconds to samples
conversion_factor = sampling_frequency/1000; % based on how many samples the data records in a second
bin_spikes_every_x_samples = bin_spikes_every_x_millisecond * conversion_factor;

%% Time bins every how many seconds?
spike_density = zeros(row_spike_identification, after_binning_data_size);

for n_electrodes=1:row_spike_identification
    n_samples_spike_density = 1;
    for n_samples=1:bin_spikes_every_x_samples:column_spike_identification
        
        if n_samples+bin_spikes_every_x_samples > column_spike_identification
            break
        
        else
            electrode_section_data = electrode_spike_identification(n_electrodes, n_samples:n_samples+bin_spikes_every_x_samples);
            find_where_spikes = find(electrode_section_data==1);
            spike_density(n_electrodes, n_samples_spike_density) = length(find_where_spikes);
        end
    n_samples_spike_density = n_samples_spike_density+1;    
    end
    fprintf("Finished electrode %d spike count.\n", n_electrodes)
end

end

