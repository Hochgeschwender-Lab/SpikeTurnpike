clc; clear;

%% Define directory to access data from particular folders
% Overarching downsampled data directory
projectFolder = uigetdir('', 'Select the project folder in which you want to analyze multiple groups');
projectFolder = fullfile(projectFolder,'SpikeStuff');

dinfo = dir(projectFolder);
dinfo(~[dinfo.isdir]) = [];  %remove non-directories
dinfo(ismember({dinfo.name}, {'.', '..'})) = [];  %remove . and .. directories
foldernames = fullfile(projectFolder, {dinfo.name});
numfolders = length(foldernames);

for ii = 1:numfolders
    [~,groupName,~] = fileparts(foldernames(ii));
    fprintf('Loading data for %s group (%d/%d)\n',groupName,ii,numfolders);

    reordered_data_directory = foldernames{ii};
    reordered_data_directory = strcat(reordered_data_directory, '/MUA/');
    %get_data_for_folders_starting_with = 'exp'; % CHANGED BY GRANT: original was 'emx'
    
    all_files_in_directory = dir(reordered_data_directory);
    logical_vector_to_classify_directories_within = [all_files_in_directory.isdir];
    get_directories_within = all_files_in_directory(logical_vector_to_classify_directories_within);
    all_folders_within_directory = {get_directories_within(3:end).name};
    folders_to_extract_data_from = {};
    
    % Initialize how many files you want to extract data from
    total_number_of_files_to_extract_from = 0;
    
    for n_folders=1:length(all_folders_within_directory)
        folder_name = all_folders_within_directory(n_folders);
        
        % Convert the folder_name from cell type to array type
        folder_name = cell2mat(folder_name);
        first_three_letters_of_folder_name = folder_name(1:3);
        
        % Not case sensitive
        %if strcmpi(first_three_letters_of_folder_name, get_data_for_folders_starting_with)
            folders_to_extract_data_from{end+1} = folder_name;
            total_number_of_files_to_extract_from = total_number_of_files_to_extract_from+1;
        %end
        
    end
    
    fprintf("This script will extract data from %d folders.\n", total_number_of_files_to_extract_from)
    % Run a loop through the available downsampled data in the specified folder
    for n_folders=1:total_number_of_files_to_extract_from
        fprintf("Working on file number %d...\n", n_folders)
        file_directory_to_save_data = strcat(reordered_data_directory, folders_to_extract_data_from(n_folders), '/');
        file_directory_to_save_data = char(file_directory_to_save_data);
        
        %{
            A folder within the specified directory should already have
            sub-folders named "allData" and "figures" from the downsampled
            section already set up. But, if not then the following code will
            make the directory. What this function basically does is that it
            checks if a folder within the "folder_directory_to_save" exists or
            not. If it does exist, then does nothing, but if the folder does
            not exist, it makes the particular folder within the specified
            directory.
        
        %}
        
        % Save all the data to a separate folder called "allData".
        check_if_folder_exists(strcat(file_directory_to_save_data, 'allData'));
        
        % Save all the figures to a separate folder called "figures".
        check_if_folder_exists(strcat(file_directory_to_save_data,  'figures'));

        
        %% Get the electrode reordering data
        electrodes_order_data_file_name = 'electrodes_order.mat';
        electrodes_order_data_directory = strcat(file_directory_to_save_data, 'allData', '/', electrodes_order_data_file_name);
        
        % Access the data by using the above defined path directory
        electrodes_order_data_struct = load(electrodes_order_data_directory);
        
        % The channel data is organized as a struct data type so 
        % the data will need to be extracted through the fields within channel
        electrodes_order = electrodes_order_data_struct.electrodes_order;

        %% Access the time stamps in microseconds
        time_stamps_file_name = 'current_timestamp_ms.mat';
        time_stamps_directory = strcat(file_directory_to_save_data, 'allData', '/', time_stamps_file_name);
    
        % Access the data by using the above defined path directory
        time_stamps_struct = load(time_stamps_directory);
    
        % The channel data is organized as a struct data type so 
        % the data will need to be extracted through the fields within channel
        time_stamps = time_stamps_struct.current_timestamps_ms; %int64
        time_stamps = double(time_stamps);

        [row_time_stamps, column_time_stamps] = size(time_stamps);
        total_trials = row_time_stamps; %this assumes 5 stims, see comment 

        % Get the number of rows and columns of the channel_data for MEA,
        % should be 5 stims using Macy's protocol, but need to treat as one
        % as it covers roughly 40ms of time. Therefore, we will pull out
        % the first onset and the last offset and treat like one trial.  

        %onset_current = time_stamps(1,1);
        %offset_current = time_stamps(end, 2);
        %time_stamps_mea = [onset_current,offset_current,1];
        %[row_time_stamps, column_time_stamps] = size(time_stamps_mea);


        
        %% Access the spike density data in us
        spike_density_data_file_name = 'spike_density_for_analyzing_change_us.mat';
        particular_spike_density_data_directory = strcat(file_directory_to_save_data, 'allData', '/', spike_density_data_file_name);
    
        % Access the data by using the above defined path directory
        spike_density_struct = load(particular_spike_density_data_directory);
    
        % The channel data is organized as a struct data type so 
        % the data will need to be extracted through the fields within channel
        spike_density = spike_density_struct.spike_density_for_analyzing_change_us;

        % Get the number of rows and columns of the channel_data
        [row_spike_density, column_spike_density] = size(spike_density);
        total_electrodes = row_spike_density;


        
        
        %% NOT NORMALIZED "spike_density" in microseconds
        % Extract NOT NORMALIZED time stamp values from "spike_density"
        before_timestamp_range = 2; % 
        after_timestamp_range = 7;% 

        extracted_time_stamps_NOT_normalized_us = ...
        extract_timestamp_values_us_fun(spike_density, before_timestamp_range, ...
            after_timestamp_range, time_stamps);

        % Save the "extracted_time_stamps_NOT_normalized" value as a .mat file.
        save(strcat(file_directory_to_save_data, 'allData/extracted_time_stamps_NOT_normalized_us.mat'), ...
        'extracted_time_stamps_NOT_normalized_us', '-v7.3')
        



    end

end

%% THE END