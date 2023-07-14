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

    downsampled_data_directory = foldernames{ii};
    downsampled_data_directory = strcat(downsampled_data_directory, '/MUA/');
    %get_data_for_folders_starting_with = 'exp'; % CHANGED BY GRANT: original is 'emx'
    
    all_files_in_directory = dir(downsampled_data_directory);
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
        %first_three_letters_of_folder_name = folder_name(1:3);
        
        % Not case sensitive
        %if strcmpi(first_three_letters_of_folder_name, get_data_for_folders_starting_with)
            folders_to_extract_data_from{end+1} = folder_name;
            total_number_of_files_to_extract_from = total_number_of_files_to_extract_from+1;
        %end
        
    end
    
    fprintf("This script will extract data from %d folders.\n", total_number_of_files_to_extract_from)
    
    % Run a loop through the available downsampled data in the specified folder
    for n_folders=1:total_number_of_files_to_extract_from
        fprintf("  Working on folder number %d...\n", n_folders)
        file_directory_to_save_data = strcat(downsampled_data_directory, folders_to_extract_data_from(n_folders), '/');
        file_directory_to_save_data = char(file_directory_to_save_data);

        if isfile(fullfile(file_directory_to_save_data,'allData',"reordered_electrode_data_downsampled.mat"))
            fprintf("    Reordered downsampled data file already exists, skipping...\n");
        else
            
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
        
            %% Access the downsampled data
            downsampled_data_file_name = 'electrode_data_downsampled.mat';
            particular_downsampled_data_directory = strcat(file_directory_to_save_data, 'allData', '/', downsampled_data_file_name);
            
            % Access the data by using the above defined path directory
            electrode_data_downsampled_struct = load(particular_downsampled_data_directory);
            
            % The channel data is organized as a struct data type so 
            % the data will need to be extracted through the fields within channel
            electrode_data_downsampled = electrode_data_downsampled_struct.electrode_data_downsampled;
            electrode_data_downsampled = electrode_data_downsampled(1:32, :); 
            %the above pulls out only the ephys electrodes 1-32, not the analog signals on electrdoes 33-35
    
            % Get the number of rows and columns of the channel_data
            [row_downsampled, column_downsampled] = size(electrode_data_downsampled);
            total_electrodes = row_downsampled;
        
            % Divide the data by 4 to maintain uniformity
            electrode_data_downsampled = electrode_data_downsampled/4;
            %% Noise removal Section
        
            %{
                What this section achives is that if the electrode rms values do not
                align within the specified criteria, then it classifies the electrodes
                as noise.
            %}
        
            fprintf("    Checking which electrodes are noisy and which are not...\n")
            [electrode_number_accept, electrode_number_reject] = noise_removal_fun(electrode_data_downsampled, file_directory_to_save_data);
            fprintf("    Finished checking which electrodes are noisy!\n")
        
            %% Reorder the electrodes as per the given order.
            %% Access the electrodes order data
            electrodes_order_data_file_name = 'electrodes_order.mat';
            electrodes_order_data_directory = strcat(file_directory_to_save_data, 'allData/', electrodes_order_data_file_name);
            
            % Access the data by using the above defined path directory
            electrodes_order_data_struct = load(electrodes_order_data_directory);
            
            % The channel data is organized as a struct data type so 
            % the data will need to be extracted through the fields within channel
            electrodes_order = electrodes_order_data_struct.electrodes_order;
            
            % Define a matrix to store the samples in the defined new order
            reordered_electrode_data_downsampled = zeros(row_downsampled, column_downsampled);
            
            fprintf("    Reordering electrodes/channels...\n")
            for n_electrodes=1:total_electrodes
                original_electrode_number = n_electrodes;
                new_electrode_number = get_electrodes_order_fun(electrodes_order, original_electrode_number, file_directory_to_save_data);
                
                reordered_electrode_data_downsampled(new_electrode_number, :) = electrode_data_downsampled(original_electrode_number, :);
                fprintf("      Electrode %d --> Electrode %d\n", original_electrode_number, new_electrode_number)
            end
            fprintf("    Finished reordering electrodes!\n")
            
            fprintf("    Saving the reordered data to file...\n")
            
            % Save the downsampled .mat data.
            save(strcat(file_directory_to_save_data, 'allData/reordered_electrode_data_downsampled.mat'),...
                'reordered_electrode_data_downsampled', '-v7.3')
            fprintf("    Finished saving the reordered data to file!\n")
        end
    end
end

%% THE END