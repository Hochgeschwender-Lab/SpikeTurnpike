function [mean_from_extracted_timestamp_data, SEM_from_extracted_timestamp_data] = ...
    SEM_calculate_fun(extracted_time_stamp_data, total_time_range, key_word, file_directory)


[row_matrix, column_matrix, depth_matrix] = size(extracted_time_stamp_data);

% Take the mean of the extracted time stamp data passed onto this function 
mean_from_extracted_timestamp_data = zeros(depth_matrix, total_time_range);
SEM_from_extracted_timestamp_data = zeros(depth_matrix, total_time_range);

for n_electrodes=1:depth_matrix
    for time_stamp=1:total_time_range
        mean_from_extracted_timestamp_data(n_electrodes, time_stamp) = ...
            mean(extracted_time_stamp_data(:, time_stamp, n_electrodes));
        SEM_from_extracted_timestamp_data(n_electrodes, time_stamp) = ...
            std(extracted_time_stamp_data(:, time_stamp, n_electrodes))/sqrt(total_time_range);
    end
end

% Save the "SEM_from_extracted_timestamp_data" value as a .mat file.
save(strcat(file_directory, 'allData/SEM_from_extracted_timestamp_', key_word, '_data.mat'), ...
    'SEM_from_extracted_timestamp_data', '-v7.3')

end

