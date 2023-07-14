clear;
turn_warnings_off = warning ('off', 'all');

%% Access the .nev file from Blackrock

projectFolder = uigetdir('', "Select a project folder where you have stored several subfolders of NEV files.");
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

%% Iterate through all the individual folders to find "nev" files
for n_folders = 1:numfolders
    
    % Find the particular file directory where the raw data ".nev" data are stored
    raw_data_directory = folder_names{n_folders};
    
    % Check if 'MUA' folder exists, if not create one!
    check_if_folder_exists(fullfile(raw_data_directory, "MUA"));
    
    [~, folder_name, ~] = fileparts(raw_data_directory);
    fprintf("Loading data for %s group (%d/%d)\n", folder_name, n_folders, numfolders);
    
    %% Get all files with extension ".ns6" and ".nev" files in specified file directory
    all_ns6_file_name_struct = dir(fullfile(raw_data_directory, "*.ns6"));
    all_nev_file_name_struct = dir(fullfile(raw_data_directory, "*.nev"));
    
    % Find the number of ".ns6" and ".nev" files
    number_of_ns6_files = length(all_ns6_file_name_struct);
    number_of_nev_files = length(all_nev_file_name_struct);
    
    fprintf("There are %d .ns6 and %d .nev files in the file directory specified.\n\n", number_of_ns6_files, number_of_nev_files);
    
    %% Start a loop to iteratie through all recognized ".ns6" and ".nev" files
    % Start a loop to iterate through all recognized ".nsx" files
    
    for file_number=1:number_of_ns6_files
        % Get the file name in every iteration file name
        ns6_file_name = all_ns6_file_name_struct(file_number).name;
        nev_file_name = all_nev_file_name_struct(file_number).name;
        
        % Define the full path of the .nev file
        nev_file_path = strcat(raw_data_directory, '/', nev_file_name);
        
        % Define a clear path where you want to save the downsampled data
        file_directory_to_save_data = fullfile(raw_data_directory, "MUA", ns6_file_name(1:end-4));
        
        fprintf("Working on file number %d titled '%s'...\n", file_number, nev_file_name)
        
        nev_file = openNEV(nev_file_path, 'read').Data;
        
        nev_serialIO_TimeStampSec = nev_file.SerialDigitalIO.TimeStampSec;
        timestamp_ms = nev_serialIO_TimeStampSec*1000; % ms
        
        check_if_folder_exists(fullfile(file_directory_to_save_data, "allData"));
        fprintf("Saving in %s\n", fullfile(file_directory_to_save_data, 'allData', 'timestamp_ms_nev.mat'));
        save(fullfile(file_directory_to_save_data, 'allData', 'timestamp_ms_nev.mat'), 'timestamp_ms');
        
    end
end

fprintf("Finished converting all .nev files to a .mat file!\n")

%% THE END