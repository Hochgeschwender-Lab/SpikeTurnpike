clear;

%% Access data from the saved timestamps from the "a2_timestamp_extraction_whisker_stim.ipynb" file

projectRootFolder = uigetdir('', "Select which project folder you are extracting timestamps from");
projectFolder = fullfile(projectRootFolder,'SpikeStuff');

% Get the directory information from the specified project folder
dinfo = dir(projectFolder);

% Remove all path names that are not directories
dinfo(~[dinfo.isdir]) = [];

% Remove '.' and ".." directories
dinfo(ismember({dinfo.name}, {'.'})) = [];
dinfo(ismember({dinfo.name}, {'..'})) = [];  

% Using the extracted directory information, collect the folder names of folders that have the ".nsx" file stored
folder_names = fullfile(projectFolder, {dinfo.name});

% Find how many folders you need to extract data from
numfolders = length(folder_names);

%% Iterate through all the individual folders to find the "without_stimulation_intensity_timestamp_ms.mat" files
for n_folders=1:numfolders
    
    % Find the particular file directory where the raw data ".nsx" data are stored
    raw_data_directory = folder_names{n_folders};
    
    % Check if 'MUA' folder exists, if not create one!
    check_if_folder_exists(fullfile(raw_data_directory, "MUA"))
    
    [~, folder_name, ~] = fileparts(raw_data_directory);
    fprintf("Loading data for %s group (%d/%d)\n", folder_name, n_folders, numfolders);
    
    %% Get all files with extension ".ns6" and ".nev" files in specified file directory
    all_ns6_file_name_struct = dir(fullfile(raw_data_directory, "*.ns6"));
    
    % Find the number of ".ns6" and ".nev" files
    number_of_ns6_files = length(all_ns6_file_name_struct);
    
    %% Start a loop to iterate through all recognized ".ns6" and ".nev" files
    for file_number=1:number_of_ns6_files
        
        % Get the file name in every iteration file name
        ns6_file_name = all_ns6_file_name_struct(file_number).name;
        
        % Define the full path of the .ns6 file
        ns6_file_path = fullfile(raw_data_directory, ns6_file_name);
        
        % Get the ns6 file to count frequencies in the analog signal
        ns6_data = openNSx('read', ns6_file_path).Data;
        analog_signal_optostim = ns6_data(33, :);
        clear ns6_data
        
        all_data_directory = fullfile(raw_data_directory, "MUA", ns6_file_name(1:end-4), "allData");
        check_if_folder_exists(all_data_directory);
        
        fprintf("Working with folder %s...", ns6_file_name(1:end-4))
        % timestamps in ms
        timestamps_without_frequency = load(fullfile(all_data_directory, "without_frequency_timestamp_ms.mat")).data;
        
        num_of_trials_to_be_removed = sum(timestamps_without_frequency(:, 1)<=0);
        timestamps_without_frequency(end-num_of_trials_to_be_removed+1:end, :) = [];
        num_of_trials = length(timestamps_without_frequency);
        
        %% Apply filter to the analog signal to count frequencies
        frequency_count_identification_smoothening_method = 'gaussian';
        
        % Store all recorded frequencies in an array
        all_frequencies = [];
        
        % specify the timewindow you want to analyze before and after stimulation
        before_stimulation = 500;
        after_stimulation = 500;
        
        % define sampling factor to convert the milliseconds to samples
        sampling_factor = 30; % conversion from ms to samples
        
        for n_trials=1:num_of_trials
            onset = timestamps_without_frequency(n_trials, 1) * sampling_factor;
            offset = timestamps_without_frequency(n_trials, 2) * sampling_factor;
            
            try
                all_frequencies(end+1) = find_number_of_frequencies_fun(onset, offset, before_stimulation, after_stimulation, ...
                    sampling_factor, frequency_count_identification_smoothening_method, analog_signal_optostim);
            catch
                timestamps_without_frequency(n_trials, :) = [];
                continue
            end

        end
        
        all_available_frequencies = unique(all_frequencies)
        
        % Save the timestamps in milliseconds
        timestamp_ms = timestamps_without_frequency;
        timestamp_ms(:, 3) = all_frequencies();

        file_directory_to_save_data = all_data_directory;
        save(fullfile(file_directory_to_save_data, 'timestamp_ms.mat'), 'timestamp_ms')
        fprintf("Saved the timestamps in milliseconds.\n")
        
        % Save the timestamps in seconds
        timestamp_s = zeros(size(timestamp_ms)); % s
        timestamp_s(:, 1:2) = timestamp_ms(:, 1:2)/1000; % s
        timestamp_s(:, 3) = timestamp_ms(:, 3); % s
        save(fullfile(file_directory_to_save_data, 'timestamp_s.mat'), 'timestamp_s')
        fprintf("Saved the timestamps in seconds.\n")
  
    
    end
    
end









