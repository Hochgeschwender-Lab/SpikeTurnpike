function downsample_data_fun(downsampling_factor, ns6_file_directory, file_directory_to_save_data)

%{
    In this section, we will reduce the sampling frequency from 30kHz to 10
    kHz. Start a while loop to only record the median value from every
    group of 3 events/samples.
%}

%% Receive data from the ns6_to_mat function.
%{
    This step calls a function "ns6_to_mat_fun" to access the ns6 data. The
    only input argument to this particular function would be the
    "file_directory" to specify where you want to save the data.
%}
fprintf("Reading file...\n")
electrode_data_original = ns6_to_mat_fun(ns6_file_directory);
    
%{
    Get the row and column information of the initial recorded data so that
    you can loop through the electrodes/channels (as rows) and the samples
    (as columns).
%}
    
[row_original, column_original] = size(electrode_data_original);
fprintf("\nThe size of the original data is %d x %d.\n", ...
    row_original, column_original)

% Create an empty matrix to store the values of the downsampled data
electrode_data_downsampled = zeros(row_original, ...
    ceil(column_original/downsampling_factor));

fprintf("Downsampling data...\n")
% Start a loop that loops through all electrodes
for n_electrodes=1:row_original
    
    %{
        Start a loop that loops through through all samples of a electrode
        skipping the "downsampled_factor" worth of next n_samples
    %}
    
    % Set the counter for n_samples_downsampled to update it accordingly
    % Might be an inefficient way to do this but I will update it later
    n_samples_downsampled = 1;
    
    for n_samples=1:downsampling_factor:column_original
        
        %{
            Set a condition where if the index threshold is exceeded, to
            break the loop. For example, in this example, if the threhold
            reaches an index that exceeds the n_samples in that particular
            electrode, it will break the loop else it runs the next block.
        
        %}
        
        if n_samples+downsampling_factor > column_original
            break
        
        %{
            Assign the median values of "downsampled_factor" set, in this
            case 3, to the zero matrix initiated for this purpose.
            ! Be careful that the downsampled matrix has 3 times lesser
            columns compared to the original set of matrix.
        %}
        
        else
            
            three_samples = electrode_data_original ...
                (n_electrodes, n_samples:n_samples+(downsampling_factor-1));
            median_three_samples = median(three_samples);
            electrode_data_downsampled(n_electrodes, n_samples_downsampled) ...
                = median_three_samples;
         
        end
        
        % Update the n_sample index for the downsampled matrix 
        n_samples_downsampled = n_samples_downsampled + 1; 
    end
    
    fprintf("Finished DownSampling for electrode %d.\n", n_electrodes)
end

% Get the size of downsampled channel_data.
[row_downsampled, column_downsampled] = size(electrode_data_downsampled);
fprintf("The size of the downsampled data is %d x %d.\n", ...
   row_downsampled, column_downsampled)
    
% Save the downsampled .mat data.
save(strcat(file_directory_to_save_data, '/allData/electrode_data_downsampled.mat'),...
    'electrode_data_downsampled', '-v7.3')

end

