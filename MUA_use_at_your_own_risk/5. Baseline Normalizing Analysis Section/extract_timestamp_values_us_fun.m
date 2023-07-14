function extracted_time_stamps = extract_timestamp_values_us_fun(spike_density, before_timestamp_range, after_timestamp_range, time_stamps);

%{
    Extract the time stamps in microseconds in binary and convert to us time stamps per
    trial
%}

% Get the row and column information for the "time_stamps" variable
[row_time_stamps, ~] = size(time_stamps); %should be us 

% Get the row and column information for the "spike_density" variable
[row_spike_density, column_spike_density] = size(spike_density); 

stimulation_exposure_duration = time_stamps(1, 2)-time_stamps(1, 1);
total_time_range = abs(before_timestamp_range)+stimulation_exposure_duration+after_timestamp_range;

% Set up a 3D array to identify particular electrode, time for analysis, and samples
extracted_time_stamps = zeros(row_time_stamps, total_time_range, row_spike_density);

for n_electrodes=1:row_spike_density
    
    for n_trials=1:row_time_stamps
        time_of_stimulation = time_stamps(n_trials, 1);
        lower_time_range = -before_timestamp_range+time_of_stimulation;
        upper_time_range = time_of_stimulation+after_timestamp_range+stimulation_exposure_duration;

            if upper_time_range >= column_spike_density
                upper_time_range = column_spike_density;
         
                total_time_range = upper_time_range-lower_time_range;
                particular_stimulation_spikes = spike_density(n_electrodes, ...
                        lower_time_range:upper_time_range);
                
                for n_samples=1:total_time_range
                    extracted_time_stamps(n_trials, n_samples, n_electrodes) = ...
                        particular_stimulation_spikes(n_samples);

                end

             else
                
                total_time_range = upper_time_range-lower_time_range;
                particular_stimulation_spikes = spike_density(n_electrodes, ...
                        lower_time_range:upper_time_range);
               
                for n_samples=1:total_time_range
                    extracted_time_stamps(n_trials, n_samples, n_electrodes) = ...
                        particular_stimulation_spikes(n_samples);
                end
                
            end
     
    end
end

end
