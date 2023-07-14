clc; clear;

% Starting with the following project folder structure...
%
% project_folder/
% ├─ SpikeStuff/
% │  ├─ Group_1/
% │  │  ├─ Recording_1/
% │  │  │  ├─ Recording1.ns6
% │  │  │  ├─ Recording1.nev
% │  │  │  ├─ Recording1_timestamps.txt
% │  │  │  ├─ Recording1.ccf
% │  │  ├─ Recording_2/
% │  │  ├─ Recording_N/
% │  ├─ Group_2/
% ├─ LFP (Optional)/
%
% For each recording folder, this script will


%% Get user input for project folder
initial_start_directory = '/Users/cresp1el/Documents/MATLAB'; 
projectFolder = uigetdir(initial_start_directory, 'Select the project folder in which you want to analyze multiple groups');
projectFolder = fullfile(projectFolder,'SpikeStuff');

% Get the directory information from the specified project folder
dinfo = dir(projectFolder);

% Remove all path names that are not directories
dinfo(~[dinfo.isdir]) = [];

% Remove '.' and '..' directories
dinfo(ismember({dinfo.name}, {'.', '..'})) = [];

% Using the extracted directory information, collect the group folder names
groupfoldernames = fullfile(projectFolder, {dinfo.name});

% Find how many folders you need to extract data from
numGroups = length(groupfoldernames);

%% Iterate through groups and recordings
for ii = 1:numGroups % enter for loop based on group 
    [~,this_group] = fileparts(groupfoldernames(ii)); %define the name of group, the tilde character, ~ , in the brackets allows you to surpress output of the filepart fx  
    fprintf('Loading data for %s group\n',this_group);

    groupDir = groupfoldernames{ii}; %define the directory of group
    dinfo2 = dir(groupDir); %define the contents of the directory
    dinfo2(~[dinfo2.isdir]) = []; 
    dinfo2(ismember({dinfo2.name},{'.','..'})) = []; %remove non directories from struct dinfo2
    recfoldernames = fullfile(groupDir, {dinfo2.name}); %create a 1xN cell array for each recording within group, N is number of recordings
    numRecordings = length(recfoldernames);   

    for jj = 1:numRecordings %loop through the recordings within a group folder
        [~,this_recording] = fileparts(recfoldernames(jj));
        fprintf('    Loading data for recording %s\n',this_recording);

        recDir = recfoldernames{jj};

        %% If MUA folder and sub-folders within MUA does not already exist, create it and proceed
        MUA_Directory = fullfile(recDir,"MUA"); %path of the MUA directory 
        
        MUA_allData_Directory = fullfile(recDir,"MUA/allData/"); %path of the allData directory
        MUA_figures_Directory = fullfile(recDir,"MUA/figures/"); %path of the figures directory 
   
        
            %% Access the downsampled data
            downsampled_data_file_name = 'electrode_data_downsampled.mat';
            particular_downsampled_data_directory = strcat(MUA_allData_Directory, downsampled_data_file_name);
            
            % Access the data by using the above defined path directory
            electrode_data_downsampled_struct = load(particular_downsampled_data_directory);
            
            % The channel data is organized as a struct data type so 
            % the data will need to be extracted through the fields within channel
            electrode_data_downsampled = electrode_data_downsampled_struct.electrode_data_downsampled;

            clear electrode_data_downsampled_struct

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
            [electrode_number_accept, electrode_number_reject] = noise_removal_fun(electrode_data_downsampled, MUA_allData_Directory);
            fprintf("    Finished checking which electrodes are noisy!\n")

            % Save the downsampled .mat data.
            save(strcat(MUA_Directory, '/allData/electrode_data_downsampled.mat'),...
                'electrode_data_downsampled', '-v7.3')
            fprintf("    Finished saving the reordered data to file!\n")

            clear electrode_data_downsampled 
       
    end
end
