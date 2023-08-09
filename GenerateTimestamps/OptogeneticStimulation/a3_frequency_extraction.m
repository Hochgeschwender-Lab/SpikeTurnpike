clear;
turn_warnings_off = warning ('off', 'all');

%% Access data from the saved timestamps from the "a2_b" file
%% Get user input for project folder
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
    
    %% Iterate through all the individual folders to find the "without_stimulation_intensity_timestamp_ms.mat" files

    for jj = 1:numRecordings %loop through the recordings within a group folder
        [~,this_recording] = fileparts(recfoldernames(jj));
        fprintf('    Loading data for %s recording (%d/%d)\n', this_recording, jj, numRecordings);

        recDir = recfoldernames{jj};
        %% Define MUA folder and sub-folders within MUA 
        MUA_Directory = fullfile(recDir,'MUA'); %path of the MUA directory         
        MUA_allData_Directory = fullfile(recDir,'MUA','allData'); %path of the allData directory
        MUA_figures_Directory = fullfile(recDir,'MUA','figures'); %path of the figures directory 
        

        ns6_file_name = dir(fullfile(recDir, "*.ns6")).name; %grab the name of the .ns6 file from the struct

        % Define the full path of the .ns6 file
        ns6_file_path = fullfile(recDir, ns6_file_name);

        % Get the ns6 file to count frequencies in the analog signal
        ns6_data = openNSx('read', ns6_file_path).Data;
       
        analog_signal_optostim = ns6_data(33, :);

        clear ns6_data

        fprintf("Working with folder %s...", ns6_file_name(1:end-4))

        % timestamps in ms
        timestamps_without_frequency = load(fullfile(MUA_allData_Directory, "without_frequency_timestamp_ms.mat")).data;

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

        all_available_frequencies = unique(all_frequencies);

        % Save the timestamps in milliseconds
        timestamp_ms = timestamps_without_frequency;
        timestamp_ms(:, 3) = all_frequencies();

        file_directory_to_save_data = MUA_allData_Directory;
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



