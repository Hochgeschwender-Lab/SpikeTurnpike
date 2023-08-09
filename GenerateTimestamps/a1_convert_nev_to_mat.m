clear;
turn_warnings_off = warning ('off', 'all');

%% Get user input for project folder
%initial_start_directory = '/home/cresp1el-local/Documents/MATLAB';
projectFolder = uigetdir('Select the project folder in which you want to analyze multiple groups');
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
    fprintf('Loading data for %s group (%d/%d)\n', this_group, ii, numGroups);


    groupDir = groupfoldernames{ii};
    dinfo2 = dir(groupDir);
    dinfo2(~[dinfo2.isdir]) = []; 
    dinfo2(ismember({dinfo2.name},{'.','..'})) = []; %remove non directories from struct dinfo2
    recfoldernames = fullfile(groupDir, {dinfo2.name}); %create a 1xN cell array for each recording within group, N is number of recordings
    numRecordings = length(recfoldernames);   

    %% Iterate through all the individual recordings to find "nev" files
    for jj = 1:numRecordings %loop through the recordings within a group folder
        [~,this_recording] = fileparts(recfoldernames(jj));
        fprintf('    Loading data for %s recording (%d/%d)\n', this_recording, jj, numRecordings);

        recDir = recfoldernames{jj};
        
        %% Define MUA folder and sub-folders within MUA 
        MUA_Directory = fullfile(recDir,'MUA'); %path of the MUA directory         
        MUA_allData_Directory = fullfile(recDir,'MUA','allData'); %path of the allData directory
        MUA_figures_Directory = fullfile(recDir,'MUA','figures'); %path of the figures directory 

        %% Check to see if timestamps have been extracted already, if so skip the recording, else proceed with extraction
        if isfile(fullfile(MUA_Directory,'allData','timestamp_ms_nev.mat'))
        
            fprintf('Timestamps have been extracted, skipping...\n')
        else
            %make the directory if it does not exsist

            check_if_folder_exists(MUA_allData_Directory); % make the MUA/allData
            check_if_folder_exists(MUA_figures_Directory); % make the MUA/figures


            %% Get file with extension ".ns6" and ".nev"  in specified directory of the recording         
            ns6_file_name = dir(fullfile(recDir, '*.ns6'));
            nev_file_name = dir(fullfile(recDir, '*.nev'));
    
            % Find the number of ".ns6" and ".nev" files
            number_of_ns6_files = length(ns6_file_name);
            number_of_nev_files = length(nev_file_name);
    
            fprintf('    There are %d .ns6 and %d .nev files in the recording directory specified.\n\n', number_of_ns6_files, number_of_nev_files)
    
            % Open NEV file, pullout timestamps in s, convert into ms, and save
            nev_file_path = fullfile(recDir, nev_file_name.name);
            nev_file = openNEV(nev_file_path, 'read').Data;
    
            nev_serialIO_TimeStampSec = nev_file.SerialDigitalIO.TimeStampSec;
            timestamp_ms = nev_serialIO_TimeStampSec*1000; % ms
    
            fprintf("Saving in %s\n", fullfile(MUA_Directory, 'allData', 'timestamp_ms_nev.mat'));
            save(fullfile(MUA_Directory, 'allData', 'timestamp_ms_nev.mat'), 'timestamp_ms')

        end

    end

end

fprintf('Finished converting all .nev files to a .mat file!\n')


