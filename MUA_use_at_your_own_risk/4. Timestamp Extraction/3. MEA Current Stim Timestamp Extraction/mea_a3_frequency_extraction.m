clc; clear;

%% Access data from the saved timestamps from the "a2_timestamp_extraction_whisker_stim.ipynb" file

projectFolder = uigetdir('', "Select which project folder you are extracting timestamps from");
projectFolder = fullfile(projectFolder,'SpikeStuff');

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
    check_if_folder_exists(strcat(raw_data_directory, "/MUA/"))
    
    [~, folder_name, ~] = fileparts(raw_data_directory)
    fprintf("Loading data for %s group (%d/%d)\n", folder_name, n_folders, numfolders);
    
    %% Get all files with extension ".h5" and ".nev" files in specified file directory
    all_h5_file_name_struct = dir(fullfile(raw_data_directory, "*.h5"));
    
    % Find the number of ".ns6" and ".nev" files
    number_of_h5_files = length(all_h5_file_name_struct);
    
    %% Start a loop to iterate through all recognized ".h5"
    for file_number=1:number_of_h5_files
        
        % Get the file name in every iteration file name
        h5_file_name = all_h5_file_name_struct(file_number).name;
        
        % Define the full path of the .ns6 file
        h5_file_path = strcat(raw_data_directory, '/', h5_file_name);
        
        % Get the h5 file to extract timestamps
        cfg = [];
        cfg.dataType = 'raw';
        h5_data = McsHDF5.McsData(h5_file_path, cfg);

        h5_data = McsHDF5.McsData('/home/cresp1el-local/Documents/MATLAB/MEA_data_pipeline/TestProject/SpikeStuff/Group_1/2022-08-15T17-23-14McsRecording_2.h5', cfg);
        current_timestamps_us = h5_data.Recording{1}.EventStream{1}.Events
           
        %pull out onset & offsets of current stimulation in h5 files, convert to ms, construct
        % .mat file where columns are onset, offset

        current_timestamps_us_onsets = h5_data.Recording{1}.EventStream{1}.Events{1}(1,:);
        current_timestamps_us_offsets = h5_data.Recording{1}.EventStream{1}.Events{2}(1,:);
        current_timestamps_us_matrix = [current_timestamps_us_onsets(:) , current_timestamps_us_offsets(:)];
        
        %convert microsecond (us) matrix into ms and s timestamps
        
        current_timestamps_us = current_timestamps_us_matrix; 
        current_timestamps_s = current_timestamps_us_matrix / 1e6;
        current_timestamps_ms = current_timestamps_us_matrix / 1000;

        [stimulation_number,~] = size(current_timestamps_us_matrix) % pull out number of stims
        current_stim_ids = repelem(1,stimulation_number) % 1 indicates a current stim

        
        current_timestamps_ms = [current_timestamps_ms, current_stim_ids(:)];
        %since the stim occured within ms of each other we will use the _ms
      
        current_timestamps_s = [current_timestamps_s, current_stim_ids(:)];
        % converting to _s may lose information if the stim occur in less than 1s 

        current_timestamps_us = [current_timestamps_us, current_stim_ids(:)];
        
        
        all_data_directory = strcat(raw_data_directory, "/MUA/", h5_file_name(1:end-3), "/allData/");
        check_if_folder_exists(all_data_directory);
        
        fprintf("Working with folder %s...", h5_file_name(1:end-3))

        file_directory_to_save_data = all_data_directory;
        
        % Save the timestamps in milliseconds
        save(strcat(file_directory_to_save_data, 'current_timestamp_ms.mat'), 'current_timestamps_ms')
        fprintf("Saved the timestamps in milliseconds.\n")
        
        % Save the timestamps in seconds
        save(strcat(file_directory_to_save_data, 'current_timestamp_s.mat'), 'current_timestamps_s')
        fprintf("Saved the timestamps in seconds.\n")

                % Save the timestamps in milliseconds
        save(strcat(file_directory_to_save_data, 'current_timestamp_us.mat'), 'current_timestamps_us')
        fprintf("Saved the timestamps in microseconds.\n")        
  
    
    end
    
end









