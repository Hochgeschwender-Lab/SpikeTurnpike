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
initial_start_directory = '/home/cresp1el-local/Documents/MATLAB';
%initial_start_directory = '/Users/cresp1el/Documents/MATLAB'; 
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
for ii = 1:numGroups
    [~,this_group] = fileparts(groupfoldernames(ii));
    fprintf('Loading data for %s group\n',this_group);

    groupDir = groupfoldernames{ii};
    dinfo2 = dir(groupDir);
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
   

        if ~exist(MUA_Directory,'dir')
            fprintf('MUA directory does not exist --> creating the MUA directory for %s\n',this_recording);
            mkdir(MUA_Directory);
        else
            fprintf('WARNING: MUA directory already exists --> skipping %s\n',this_recording);
        end

        if ~exist(MUA_allData_Directory,'dir')
            fprintf('allData sub-directory does not exist --> creating the MUA directory for %s\n',this_recording);
            mkdir(MUA_allData_Directory)

            %% Open NSx file and write data to 
            NSx_file = dir(fullfile(recDir,'*.ns6'));
            NSx_filepath = fullfile(recDir, NSx_file.name); 
    
            %% Specify electrode order and save .mat file for later usage
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
            save(strcat(MUA_allData_Directory,...
            'electrodes_order.mat'), 'electrodes_order', '-v7.3');

            %% Downsampling Section.
            %{
                In this section, we will reduce the sampling frequency from 30kHz to 10
                kHz. Start a while loop to only record the median value from every
                group of 3 events/samples.
            %}
        
            % Define how many sample section we want to compute the median from.
            downsampling_factor = 3;
        
            % Run the "downsample_data_fun" to downsample the data from original data.
            
            if isfile(fullfile(MUA_Directory,'allData',"electrode_data_downsampled.mat"))
            fprintf("    Downsampled data file already exists, skipping...\n");
            else
            tic
            downsample_data_fun(downsampling_factor, NSx_filepath, MUA_Directory);
            toc
            fprintf("    Finished Data downsampling!\n")
            end

        else 
            fprintf('WARNING: allData sub-directory already exists --> skipping %s\n',this_recording);
        end

        if ~exist(MUA_figures_Directory,'dir')
            fprintf('figures sub-directory does not exist --> creating the MUA directory for %s\n',this_recording);
            mkdir(MUA_figures_Directory)
        else 
            fprintf('WARNING: figures sub-directory already exists --> skipping %s\n',this_recording);
        end
    end
end
