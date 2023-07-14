function [normalized_spike_counts, normalizing_mean_array] = normalize_spike_density_fun(spike_density, extracted_time_stamps_NOT_normalized, ...
    normalized_spike_counts, before_timestamp_range, after_timestamp_range, time_stamps, file_directory)

% Get the row and column information for the "time_stamps" variable
[row_time_stamps, column_time_stamps] = size(time_stamps);

% Get the row and column information for the "spike_density" variable
[row_spike_density, column_spike_density] = size(spike_density);

stimulation_exposure_duration = time_stamps(1, 2)-time_stamps(1, 1);

normalizing_mean_array = zeros(row_spike_density, row_time_stamps);

for n_electrodes=1:row_spike_density
    for n_trials=1:row_time_stamps
        time_of_stimulation = time_stamps(n_trials, 1);
        lower_time_range = -before_timestamp_range+time_of_stimulation;
        upper_time_range = time_of_stimulation+after_timestamp_range+stimulation_exposure_duration-1;
        
        if upper_time_range > column_spike_density
            upper_time_range = column_spike_density;
            time_span_until_stimulation = length(1:before_timestamp_range);
            % time_span = length(lower_time_range:upper_time_range);

            normalizing_mean = mean(extracted_time_stamps_NOT_normalized(n_trials, 1:time_span_until_stimulation, n_electrodes));
            normalizing_mean_array(n_electrodes, n_trials) = normalizing_mean;

            normalized_spike_counts(n_electrodes, lower_time_range:upper_time_range) = ...
                spike_density(n_electrodes, lower_time_range:upper_time_range)/normalizing_mean;
 
        else
            time_span_until_stimulation = length(1:before_timestamp_range);
            % time_span = length(lower_time_range:upper_time_range);
                
            normalizing_mean = mean(extracted_time_stamps_NOT_normalized(n_trials, 1:time_span_until_stimulation, n_electrodes));
            normalizing_mean_array(n_electrodes, n_trials) = normalizing_mean;
                
            normalized_spike_counts(n_electrodes, lower_time_range:upper_time_range) = ...
                spike_density(n_electrodes, lower_time_range:upper_time_range)/normalizing_mean;
        end
    end
end

% % Save the "spike_density" as a .mat data.
% save(strcat(file_directory, 'allData/spike_density.mat'),...
%     'spike_density', '-v7.3')

end

