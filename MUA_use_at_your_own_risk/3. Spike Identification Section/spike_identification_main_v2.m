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
   
        %% Get the electrode reordering data
        electrodes_order_data_file_name = 'electrodes_order.mat';
        electrodes_order_data_directory = strcat(MUA_allData_Directory,electrodes_order_data_file_name);
        
        % load the data by using the above defined path directory
        electrodes_order_data_struct = load(electrodes_order_data_directory);
        
        % The channel data is organized as a struct data type so 
        % the data will need to be extracted through the fields within channel
        electrodes_order = electrodes_order_data_struct.electrodes_order;
        
        %% Access the reordered downsampled data
        reordered_data_file_name = 'electrode_data_downsampled.mat';
        particular_reordered_data_directory = strcat(MUA_allData_Directory, reordered_data_file_name);
        
        % Access the data by using the above defined path directory
        electrode_data_reordered_struct = load(particular_reordered_data_directory);
        
        % The channel data is organized as a struct data type so 
        % the data will need to be extracted through the fields within channel
        electrode_data_reordered = electrode_data_reordered_struct.electrode_data_downsampled;
    
        % Get the number of rows and columns of the channel_data
        [row_reordered, column_reordered] = size(electrode_data_reordered);
        total_electrodes = row_reordered;
        
        %% Filter the data to identify spikes
        % Butterworth bandpass filter design
        fprintf("\nButterworth BandPass Filtering electrode data in progress...\n")
        % This step takes roughly 20-21 minutes to complete.
        filtered_electrode_data = zeros(size(electrode_data_reordered));
        
        % Define low pass and stop frequency
        low_frequency_pass = 500; % Hz % From Manny's code.
        low_frequency_stop = 300; % Hz % From Prakash et al.
    
        % Define high cutoff frequency
        high_frequency_stop = 4000; % Hz % From Manny's code.
        high_frequency_pass = 3000; % Hz % From Prakash et al.
    
        % Define sampling_frequency
        sampling_frequency = 10e3; % Hz % Downsampled from 30kHz
    
        % Find nyquist_frequency for normalizing the frequency
        nyquist_frequency = (sampling_frequency/2);
    
        % Determine passband and stopband edge frequency
        Wp = [low_frequency_pass high_frequency_pass] / nyquist_frequency;
        Ws = [low_frequency_stop high_frequency_stop] / nyquist_frequency;
    
        % Assign dB of ripple and attenuation
        [N, Wn] = buttord(Wp, Ws, 3, 20);
        [B, A] = butter(N, Wn);
        
        for n_electrodes=1:total_electrodes
            filtered_electrode_data(n_electrodes, :) = filtfilt(B, A, electrode_data_reordered(n_electrodes, :));
            fprintf('Filtered Electrode %d\n', n_electrodes)
        end
       
        fprintf("Finished filtering electrode data!\n")
        clear electrode_data_reordered;

        % fprintf("Saving filtered_electrode_data...\n")
        % save(strcat(MUA_Directory, '/allData/filtered_electrode_data.mat'), 'filtered_electrode_data', '-v7.3')
    
        %{
            Plot the figures separately to observe how the signals are filtered and
            set up a loop to save all the graphs.
        
        
        fprintf("Plotting data before and after filter application...\n")
        for n_electrodes=1:total_electrodes
            new_electrode_number = n_electrodes;
            original_electrode_number = electrodes_order(new_electrode_number);
            plot_filtered_data_fun(filtered_electrode_data(new_electrode_number, :), electrode_data_reordered(new_electrode_number, :), original_electrode_number, file_directory_to_save_data);
        end
        %}
        
        
        %% Spike identification section
        fprintf("\nIdentifying spikes...\n")
        electrode_spike_identification = spike_identification_fun(filtered_electrode_data, MUA_allData_Directory);
        fprintf("Finished identifying spikes!\n\n")
        
        clear filtered_electrode_data;
        %% Bin the identified spikes into 1 second intervals
    
        %{
            Binning the spike counts in a second effectively makes the 
            "spike_density" variable "firing rate".
        %}
        
        % Bin spikes for further data analysis
        % Make time bins for every 'x' second
        bin_spikes_every_x_millisecond = 1; % ms
       
        % Get the "spike_density" data from the "spike_density_fun" function
        fprintf('Binning the spikes in %d ms bins for further MUA analysis...\n', bin_spikes_every_x_millisecond)
        
        spike_density_for_analyzing_change = spike_density_fun(electrode_spike_identification, ...
            sampling_frequency, bin_spikes_every_x_millisecond, MUA_allData_Directory);
    
        % Save the "spike_density_for_analyzing_change" as a .mat data.
        fprintf("Saving spike density data for further analysis...\n")
        save(strcat(file_directory_to_save_data, 'allData/spike_density_for_analyzing_change.mat'),...
            'spike_density_for_analyzing_change', '-v7.3')
    end
end

        
%         % Bin spikes into 1 second bins for electrode observation further data analysis
%         % Make time bins for every 'x' second
%         bin_spikes_every_x_millisecond = 1000; % ms
%         
%         fprintf('\nBinning the spikes in %d ms bins for electrode observations...\n', bin_spikes_every_x_millisecond)
%     
%         % Get the "spike_density" data from the "spike_density_fun" function
%         spike_density_for_electrode_observation = spike_density_fun(electrode_spike_identification, ...
%             sampling_frequency, bin_spikes_every_x_millisecond, file_directory_to_save_data);
        
%         fprintf("\nPlotting spike density graphs for electrode observation...\n")
%     
%         % Plot the spike density as a heat map
%         plot_heatmap_fun(spike_density_for_electrode_observation, 'SpikeDensity_HeatMap', file_directory_to_save_data)
%         
%         % Plot the spike density for each electrode at the same time
%         figure('visible', 'off');
%         tiledlayout('flow');
%         for n_electrodes=1:total_electrodes
%             original_electrode_number = n_electrodes;
%             new_electrode_number = find(electrodes_order == original_electrode_number);
%             
%             nexttile;
%             plot(spike_density_for_electrode_observation(new_electrode_number, :), '-');
%             xlabel('Electrode ' + string(original_electrode_number))
%             hold on;
%         end
%         hold off;
%     
%         filename_fig = strcat(file_directory_to_save_data, 'figures/AllElectrodes/', 'SpikeDensity_Plot', '.fig');
%         saveas(gcf, filename_fig)
%         
%         filename_pdf = strcat(file_directory_to_save_data, 'figures/AllElectrodes/', 'SpikeDensity_Plot', '.pdf');
%         exportgraphics(gcf, filename_pdf)
%         
%         % Save all the figures to a separate folder called "figures/SpikeDensity".
%         check_if_folder_exists(strcat(file_directory_to_save_data, 'figures/SpikeDensity/'));
%     
%         % Plot the spike density for individual electrodes
%         for n_electrodes=1:total_electrodes
%             figure('visible', 'off');
%             
%             new_electrode_number = n_electrodes;
%             original_electrode_number = electrodes_order(new_electrode_number); 
%     
%             plot(spike_density_for_electrode_observation(new_electrode_number, :))
%             xlabel('Original Electrode ' + string(original_electrode_number))
%             ylabel('Spike Counts binned to 1 second interval')
%             filename_png = strcat(file_directory_to_save_data, 'figures/SpikeDensity/', 'Electrode', string(original_electrode_number), '_SpikeDensity.png');
%             saveas(gcf, filename_png);
%     
%             filename_fig = strcat(file_directory_to_save_data, 'figures/SpikeDensity/', 'Electrode', string(original_electrode_number), '_SpikeDensity.fig');
%             saveas(gcf, filename_fig);
%         end
%     end
% 
% end

%% THE END