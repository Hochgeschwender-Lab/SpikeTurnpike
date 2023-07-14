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
    
    %% Get all files with the extension ".h5" in the specified file_directory
    all_h5_file_name_struct = dir(fullfile(raw_data_directory, '*.h5'));
    number_of_files = length(all_h5_file_name_struct);
    
    fprintf("There are %d .h5 files in the file directory specified.\n", number_of_files)
    
    for file_number=1:number_of_files
        
        file_name = all_h5_file_name_struct(file_number).name;
        file_number_directory = strcat(all_h5_file_name_struct(file_number).folder, '/', file_name);
        file_directory_to_save_data = strcat(downsampled_data_to_be_saved_directory, '/MUA/', file_name(1:end-3), '/');
    
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
            arrangement of the MEA electrodes on the MEA culture dish. Reorder the 
            electrodes based on the order provided and picking out particular 
            electrode sample.
        %}
        

        electrodes_order = [47; 48; 46; 45; 38; 37; 28; 36; 27; 17; 26; 16; 35; 25; 15;... 
                                14; 24; 34; 13; 23; 12; 22; 33; 21; 32; 31; 44; 43; 41;... 
                                42; 52; 51; 53; 54; 61; 62; 71; 63; 72; 82; 73; 83; 64;... 
                                74; 84; 85; 75; 65; 86; 76; 87; 77; 66; 78; 67; 68; 55;... 
                                56; 58; 57];
    
        % Save the electrodes_order for future use
        save(strcat(file_directory_to_save_data,...
            'allData/electrodes_order.mat'), 'electrodes_order', '-v7.3');
        
        %% Downsampling Section.
        %{
            In this section, we will reduce the sampling frequency from 30kHz to 10
            kHz. However, this code is for MEA so it is already at 10kHz. No downsmapling is needed. 
            Start a while loop to only record the median value from every
            group of 3 events/samples for the in vivo data. 
        %}
        
        % Define how many sample section we want to compute the median
        % from. If downsampling_factor = 1, then no median is calculated
        downsampling_factor = 1;
        
        % Run the "downsample_data_fun" to downsample the data from original data.
        if isfile(fullfile(file_directory_to_save_data,'allData',"electrode_data_downsampled.mat"))
            fprintf("    Downsampled data file already exists, skipping...\n");
        else
            tic
            mea_downsample_data_fun(downsampling_factor, file_number_directory, file_directory_to_save_data);
            toc
            fprintf("    Finished Data downsampling!\n")
        end
        
    end
end

%% THE END