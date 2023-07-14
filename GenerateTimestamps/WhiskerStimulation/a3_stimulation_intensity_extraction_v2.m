clear;
turn_warnings_off = warning ('off', 'all');

%% Access data from the saved timestamps from the "a2_timestamp_extraction_whisker_stim.ipynb" file
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
        
        % timestamps in ms
        timestamps_without_stimulation_intensity = load(fullfile(MUA_Directory,'allData','without_stimulation_intensity_timestamp_ms.mat')).data;
        
        %% Access the "print_lines" file to extract stimulation intensity information from the ".txt" file
        all_stimulation_intensities_in_order = [];
        session_print_lines = load(fullfile(MUA_Directory,'allData','session_print_lines.mat')).data;
        num_of_trials = length(session_print_lines);
        
        %{
            Run a loop through every 10th array because the stimulation duration is
            10s. Run this through all the data to read the "words" that signify
            catch, low, med, and max levels of stimulation.
        %}
        %% CHANGE IF THE STIMULATION DURATION IS NOT THIS VALUE - GRANT
        every_x_interval_skip = 11; % steps
        
        for n_timestamps=1:every_x_interval_skip:num_of_trials
            
            %{
                Split the character arrays if there is a space between words. Do this for both on-sets and off-sets. The
                above operation results in a cell with 2 contents:  {on-set/off-set time} and {stimulation intensity} 
                both of which can be indexed. Using the stimulation intensity check if the stimulation is 'catch', 'low', 
                'med', or 'max'. If the condition satisifes, add the on-set/off-set values to the appropriate
                on-set/off-set variable.
            %}
            
            string_from_print_lines_split = split(session_print_lines(n_timestamps, :));
            stimulation_intensity = string_from_print_lines_split(2);
            
            if strcmpi(stimulation_intensity, 'catch')
                all_stimulation_intensities_in_order(end+1) = 1;
                
            elseif strcmpi(stimulation_intensity, 'low')
                all_stimulation_intensities_in_order(end+1) = 2;
                
            elseif strcmpi(stimulation_intensity, 'med')
                all_stimulation_intensities_in_order(end+1) = 3;
                
            elseif strcmpi(stimulation_intensity, 'max')
                all_stimulation_intensities_in_order(end+1) = 4;
                
            end
            
        end
        
        timestamps_with_stimulation_intensity = zeros(size(timestamps_without_stimulation_intensity));
        timestamps_with_stimulation_intensity(:, 1:2) = timestamps_without_stimulation_intensity(:, 1:2);
        timestamps_with_stimulation_intensity(:, 3) = all_stimulation_intensities_in_order;

        %% Count number of zeroes at the end to remove unidentifiable trials
        num_of_trials_to_be_removed = sum(timestamps_with_stimulation_intensity(:, 1)==0);
        timestamps_with_stimulation_intensity(end-num_of_trials_to_be_removed+1:end, :) = [];
        
        % Save the timestamps in milliseconds
        timestamp_ms = timestamps_with_stimulation_intensity; % ms
        save(fullfile(MUA_Directory,'allData','timestamp_ms.mat'), 'timestamp_ms')
        fprintf("Saved the timestamps in milliseconds.\n")
        
        % Save the timestamps in seconds
        timestamp_s = zeros(size(timestamp_ms)); % s
        timestamp_s(:, 1:2) = timestamp_ms(:, 1:2)/1000; % s
        timestamp_s(:, 3) = timestamp_ms(:, 3); % s
        save(fullfile(MUA_Directory,'allData','timestamp_s.mat'), 'timestamp_s')
        fprintf("Saved the timestamps in seconds.\n")
        
    end
    
end