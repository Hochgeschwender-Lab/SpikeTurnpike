clc; clear;

%% Define directory to save data and graphs

% Define where the raw data is stored
projectFolder = uigetdir('', 'Select a project folder where you have stored several subfolders of NSx files.');
projectFolder = fullfile(projectFolder,'SpikeStuff');

dinfo = dir(projectFolder);
dinfo(~[dinfo.isdir]) = [];  %remove non-directories
dinfo(ismember({dinfo.name}, {'.', '..'})) = [];  %remove . and .. directories
foldernames = fullfile(projectFolder, {dinfo.name});
numfolders = length(foldernames);

for n_folders = 1:numfolders
    raw_data_directory = foldernames{n_folders};
    downsampled_data_to_be_saved_directory = raw_data_directory;

    fprintf('Loading data for %s group (%d/%d)\n',fileparts(foldernames(n_folders)), n_folders,numfolders);
    
    %% Get all files with the extension ".ns6" in the specified file_directory
    all_ns6_file_name_struct = dir(fullfile(raw_data_directory, '*.ns6'));
    number_of_files = length(all_ns6_file_name_struct);
    
    fprintf("There are %d .ns6 files in the file directory specified.\n", number_of_files)
    
    for file_number=1:number_of_files
        
        file_name = all_ns6_file_name_struct(file_number).name;
        file_number_directory = strcat(all_ns6_file_name_struct(file_number).folder, '/', file_name);
        file_directory_to_save_data = strcat(downsampled_data_to_be_saved_directory, '/MUA/', file_name(1:end-4), '/');
    
        fprintf("  Working on file number %d titled '%s'...\n", file_number, file_name)
    
        %% Create sub-folders to store data
        %{
            Create a separate folder for the data to be saved. What this function
            basically does is that it checks if a folder within the
            "file_directory" exists or not. If it does exist, then does nothing,
            but if the folder does not exist, it makes the particular folder within
            the specified directory.
        %}
        
        % Save all the data to a separate folder called "allData".
        check_if_folder_exists(strcat(file_directory_to_save_data, 'allData'));
        
        % Save all the figures to a separate folder called "figures".
        check_if_folder_exists(strcat(file_directory_to_save_data,  'figures'));
        
        %% Specify which particular electrode you want to .
        %{
            This is the electrode number before we reorder them based on the
            arrangement of electrodes on the cylindrical probe. Reorder the 
            electrodes based on the order provided and picking out particular 
            electrode sample.
        %}
        
        electrodes_order = [14; 20; 16; 18; 1;31; 3; 29; 5; 27; 7; 25; 9; 23;...
                                11; 21; 13; 19; 15; 17; 12; 22; 10; 24; 8; 26;...
                                    6; 28; 4; 30; 2; 32];
    
        % Save the electrodes_order for future use
        save(strcat(file_directory_to_save_data,...
            'allData/electrodes_order.mat'), 'electrodes_order', '-v7.3');
        
        %% Downsampling Section.
        %{
            In this section, we will reduce the sampling frequency from 30kHz to 10
            kHz. Start a while loop to only record the median value from every
            group of 3 events/samples.
        %}
        
        % Define how many sample section we want to compute the median from.
        downsampling_factor = 3;
        
        % Run the "downsample_data_fun" to downsample the data from original data.
        if isfile(fullfile(file_directory_to_save_data,'allData',"electrode_data_downsampled.mat"))
            fprintf("    Downsampled data file already exists, skipping...\n");
        else
            tic
            downsample_data_fun(downsampling_factor, file_number_directory, file_directory_to_save_data);
            toc
            fprintf("    Finished Data downsampling!\n")
        end
        
    end
end

%% THE END