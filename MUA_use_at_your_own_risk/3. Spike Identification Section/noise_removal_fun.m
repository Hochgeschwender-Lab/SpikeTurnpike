function [electrode_number_accept, electrode_number_reject] = noise_removal_fun(electrode_data_downsampled, file_directory)

%% Get row and column information of the channel_data
[row_downsampled, column_downsampled] = size(electrode_data_downsampled);

%% Noise removal process
%{
    To filter the noise in the data, we will use the method demonstrated by
    Prakash et al (2022). We will first calculate the rms value of a
    channel. If the rms value is greater than 3 times the interquartile
    range above the 3rd quartile, or if the rms value is less than 3 times
    less than the interquartile range below the 1st quartile, the electrode
    contacts were marked as noise.
%}
        
%{
    Now, remove the rms value voltages that do not follow the criteria
    defined above and to do this define a loop to go through every sample
    of all channels.
%}

%{
    Define arrays to record the electrode (or channels) to see which one 
    classifies as noise
%}

electrode_number_accept = [];
electrode_number_reject = [];

%{
    Run a loop through the downsampled matrix channels to find rms values
    of each channel to check if the rms values fall within:
    first_quartile - 3 * inter_quartile_range < rms_value_channel <
    third_quartile + 3 * inter_quartile_range
%}

for n_channels=1:row_downsampled
    
    % Sort the downsampled data first to find the quartiles
    sorted_channel_data_downsampled = ...
        sort(electrode_data_downsampled(n_channels, :));

    % Find the median of the channel data
    median_channel_data = median(sorted_channel_data_downsampled);

    %{
        FIRST QUARTILE ALGORITHM:
        1. Find the data in the channel that falls below the median. This
        step will give you the indexes of the array values which fall below
        the median. You should expect these index to be in ascending order
        considering you just sorted the entire data in ascending order.
        2. After finding the indexes, extract the actual channel data of
        those particular indexes.
        3. From these extracted values, find the median again which will
        result in you finding the first quartile.
    %}

    % Step 1
    index_of_values_below_median = ...
        find(sorted_channel_data_downsampled < median_channel_data);

     % Step 2
    below_median_channel_data = ...
        sorted_channel_data_downsampled(index_of_values_below_median);

    % Step 3
    first_quartile = median(below_median_channel_data);

    %{
        THIRD QUARTILE ALGORITHM:
        1. Find the data in the channel that are above the median. This
        step will give you the indexes of the array values which are above
        the median. You should expect these index to be in ascending order
        considering you just sorted the entire data in ascending order.
        2. After finding the indexes, extract the actual channel data of
        those particular indexes.
        3. From these extracted values, find the median again which will
        result in you finding the third quartile.
    %}

    % Step 1
    index_of_values_above_median = ...
        find(sorted_channel_data_downsampled > median_channel_data);

    % Step 2
    above_median_channel_data = ...
        sorted_channel_data_downsampled(index_of_values_above_median);

    % Step 3
    third_quartile = median(above_median_channel_data);

    inter_quartile_range = third_quartile - first_quartile;
    
    % Find the rms value of the the particular channel
    rms_value_channel = rms(electrode_data_downsampled(n_channels, :));
    
    if (rms_value_channel > (first_quartile - 3 * inter_quartile_range)) ...
            && (rms_value_channel < (third_quartile + 3 * inter_quartile_range))
        fprintf("E%d: rms value is %0.2f < %0.2f < %0.2f: NOT NOISE\n", n_channels, ...
            first_quartile - 3 * inter_quartile_range, rms_value_channel, ...
                third_quartile + 3 * inter_quartile_range);
        electrode_number_accept = [electrode_number_accept, n_channels];
    else
        fprintf("E%d: rms value is %0.2f < %0.2f < %0.2f: NOISE\n", n_channels, ...
            first_quartile - 3 * inter_quartile_range, rms_value_channel, ...
                third_quartile + 3 * inter_quartile_range);
        electrode_number_reject = [electrode_number_reject, n_channels];
    end
end

% Save the electrodes to be accepted or rejected as .mat data.
save(strcat(file_directory, 'allData/electrode_number_accept.mat'), 'electrode_number_accept', '-v7.3')

% Save the electrodes to be rejected as .mat data.
save(strcat(file_directory, 'allData/electrode_number_reject.mat'), 'electrode_number_reject', '-v7.3')

end