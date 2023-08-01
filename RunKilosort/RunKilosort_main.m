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
% ├─ KS_config.m ***
% ├─ LFP (Optional)/
%
% *** Copy and paste KS_config.m from the folder this script is in
% (RunKilosort) to your project folder, then edit it with desired
% parameters. The point of this is reproducibility for the project.
%
% For each recording folder, this script will
%   1. Convert the .ns6 into a .bin file
%   2. Create a Kilosort/Phy working and output directory /SUA/
%   3. Run Kilosort with our standard parameters
%   4. Delete the .bin file as it is no longer needed.
%
% Then, you can use Phy to curate units in each recording.
%
% Recordings which have already been processed in Kilosort will be skipped,
% so new ones can be added without disrupting the workflow.

%% Get user input for project folder
projectRootFolder = uigetdir('Select your project folder:');
projectFolder = fullfile(projectRootFolder,'SpikeStuff');

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

%% Get number of channels
prompt = {'Enter the number of channels/electrodes on your probe:'};
dlgtitle = 'Num. Channels';
definput = {'32'};
dims = [1 40];
n_channels = inputdlg(prompt,dlgtitle,dims,definput);
n_channels = str2double(n_channels{1});
channelsArg = strcat('c:1:',num2str(n_channels)); % for openNSx call

%% Get config and channel map files
% get absolute path to Kilosort2/configFiles, where channel maps are stored
thisScipt_filepath = fileparts(mfilename('fullpath'));
idx = strfind(thisScipt_filepath,filesep);
ks_configFiles_dir = fullfile(thisScipt_filepath(1:idx(end)), 'external', 'Kilosort2', 'configFiles');
clear idx thisScipt_filepath

% prompt user to choose a specific channel map
chanMapFile = uigetfile(ks_configFiles_dir, 'Select a channel map (.mat):', fullfile(ks_configFiles_dir,'Neuronexus32channelmap_kilosortChanMap.mat'));
%chanMapFile = 'Neuronexus32channelmap_kilosortChanMap.mat';

%% Build ops struct for Kilosort
run(fullfile(projectRootFolder, 'KS_config.m'));
ops.NchanTOT  = n_channels; % total number of channels in your recording
ops.chanMap = fullfile(ks_configFiles_dir, chanMapFile);

%% Iterate through groups and recordings
for ii = 1:numGroups
    [~,this_group] = fileparts(groupfoldernames(ii));
    fprintf('Loading data for %s group\n',this_group);

    groupDir = groupfoldernames{ii};
    dinfo2 = dir(groupDir);
    dinfo2(~[dinfo2.isdir]) = [];
    dinfo2(ismember({dinfo2.name},{'.','..'})) = [];
    recfoldernames = fullfile(groupDir, {dinfo2.name});
    numRecordings = length(recfoldernames);

    for jj = 1:numRecordings
        [~,this_recording] = fileparts(recfoldernames(jj));
        fprintf('    Loading data for recording %s\n',this_recording);

        recDir = recfoldernames{jj};

        %% If SUA folder does not already exist, create it and proceed
        SUA_Directory = fullfile(recDir,"SUA");
        if ~isfile(fullfile(SUA_Directory,'params.py'))
            if ~exist(SUA_Directory,'dir')
                mkdir(SUA_Directory);
            end

            %% Open NSx file and write data to temp.bin
            NSx_file = dir(fullfile(recDir,'*.ns6'));
            NSx_filepath = fullfile(recDir, NSx_file.name);

            ns = openNSx('read', NSx_filepath, channelsArg); % changed to extract ephys channels 1:32, not 33,34,35 analog channels

            bin_filepath = fullfile(recDir,'temp.bin');
    
            if (path ~= 0)
                file = fopen(bin_filepath,'w');
                fwrite(file, ns.Data, 'int16');
                fclose(file);
            end

            clear ns

            %% Run Kilosort
            %addpath(genpath('/home/arawa1rj-local/Documents/ephGit/Kilosort2')) % path to kilosort folder
            %addpath('D:\GitHub\npy-matlab') % for converting to Phy
            rootZ = recDir; % the raw data binary file is in this folder
            rootH = SUA_Directory; % path to temporary binary file (same size as data, should be on fast SSD)
            
            ops.fproc   = fullfile(rootH, 'temp_wh.dat');
            
            % find the binary file
            %fs          = [dir(fullfile(rootZ, '*.bin')) dir(fullfile(rootZ, '*.dat'))];
            %ops.fbinary = fullfile(rootZ, fs(1).name);
            ops.fbinary = bin_filepath;
            
            rez                = preprocessDataSub(ops);
            rez                = datashift2(rez, 1);
            
            [rez1, st3, tF]     = extract_spikes(rez);
            
            rez                = template_learning(rez1, tF, st3);
            rez_save = rez1;
            tF_save = tF;
            st3_save = st3;
            
            [rez, st3, tF]     = trackAndSort(rez);
            
            rez                = final_clustering(rez, tF, st3);
            
            rez                = find_merges(rez, 1);
            
%             rootZ = fullfile(rootZ, 'kilosort3');
%             mkdir(rootZ)
            rezToPhy2(rez, SUA_Directory);

            %% Delete temp.bin to save space
            delete(fullfile(recDir,'temp.bin'));

        else
            fprintf('WARNING: SUA directory already exists --> skipping %s\n',this_recording);
        end

    end

end