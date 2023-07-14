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

    reordered_data_directory = foldernames{ii};
    reordered_data_directory = strcat(reordered_data_directory, '/MUA/');
    %get_data_for_folders_starting_with = 'exp'; % CHANGED BY GRANT: original was 'emx'
    
    all_files_in_directory = dir(reordered_data_directory);
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
        first_three_letters_of_folder_name = folder_name(1:3);
        
        % Not case sensitive
        %if strcmpi(first_three_letters_of_folder_name, get_data_for_folders_starting_with)
            folders_to_extract_data_from{end+1} = folder_name;
            total_number_of_files_to_extract_from = total_number_of_files_to_extract_from+1;
        %end
        
    end
    
    fprintf("This script will extract data from %d folders.\n", total_number_of_files_to_extract_from)
    % Run a loop through the available downsampled data in the specified folder
    for n_folders=1:total_number_of_files_to_extract_from
        fprintf("Working on file number %d...\n", n_folders)
        file_directory_to_save_data = strcat(reordered_data_directory, folders_to_extract_data_from(n_folders), '/');
        file_directory_to_save_data = char(file_directory_to_save_data);
        
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
        
        %% Get the electrode reordering data
        electrodes_order_data_file_name = 'electrodes_order.mat';
        electrodes_order_data_directory = strcat(file_directory_to_save_data, 'allData', '/', electrodes_order_data_file_name);
        
        % Access the data by using the above defined path directory
        electrodes_order_data_struct = load(electrodes_order_data_directory);
        
        % The channel data is organized as a struct data type so 
        % the data will need to be extracted through the fields within channel
        electrodes_order = electrodes_order_data_struct.electrodes_order;
        
        %% Access the reordered downsampled data
        reordered_data_file_name = 'reordered_electrode_data_downsampled.mat';
        particular_reordered_data_directory = strcat(file_directory_to_save_data, 'allData', '/', reordered_data_file_name);
        
        % Access the data by using the above defined path directory
        electrode_data_reordered_struct = load(particular_reordered_data_directory);
        
        % The channel data is organized as a struct data type so 
        % the data will need to be extracted through the fields within channel
        electrode_data_reordered = electrode_data_reordered_struct.reordered_electrode_data_downsampled;
    
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
        % save(strcat(file_directory_to_save_data, 'allData/filtered_electrode_data.mat'), 'filtered_electrode_data', '-v7.3')
    
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
        electrode_spike_identification = spike_identification_fun(filtered_electrode_data, file_directory_to_save_data);
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
            sampling_frequency, bin_spikes_every_x_millisecond, file_directory_to_save_data);
    
        % Save the "spike_density_for_analyzing_change" as a .mat data.
        fprintf("Saving spike density data for further analysis...\n")
        save(strcat(file_directory_to_save_data, 'allData/spike_density_for_analyzing_change.mat'),...
            'spike_density_for_analyzing_change', '-v7.3')
        
        % Bin spikes into 1 second bins for electrode observation further data analysis
        % Make time bins for every 'x' second
        bin_spikes_every_x_millisecond = 1000; % ms
        
        fprintf('\nBinning the spikes in %d ms bins for electrode observations...\n', bin_spikes_every_x_millisecond)
    
        % Get the "spike_density" data from the "spike_density_fun" function
        spike_density_for_electrode_observation = spike_density_fun(electrode_spike_identification, ...
            sampling_frequency, bin_spikes_every_x_millisecond, file_directory_to_save_data);
        
        fprintf("\nPlotting spike density graphs for electrode observation...\n")
    
        % Plot the spike density as a heat map
        plot_heatmap_fun(spike_density_for_electrode_observation, 'SpikeDensity_HeatMap', file_directory_to_save_data)
        
        % Plot the spike density for each electrode at the same time
        figure('visible', 'off');
        tiledlayout('flow');
        for n_electrodes=1:total_electrodes
            original_electrode_number = n_electrodes;
            new_electrode_number = find(electrodes_order == original_electrode_number);
            
            nexttile;
            plot(spike_density_for_electrode_observation(new_electrode_number, :), '-');
            xlabel('Electrode ' + string(original_electrode_number))
            hold on;
        end
        hold off;
    
        filename_fig = strcat(file_directory_to_save_data, 'figures/AllElectrodes/', 'SpikeDensity_Plot', '.fig');
        saveas(gcf, filename_fig)
        
        filename_pdf = strcat(file_directory_to_save_data, 'figures/AllElectrodes/', 'SpikeDensity_Plot', '.pdf');
        exportgraphics(gcf, filename_pdf)
        
        % Save all the figures to a separate folder called "figures/SpikeDensity".
        check_if_folder_exists(strcat(file_directory_to_save_data, 'figures/SpikeDensity/'));
    
        % Plot the spike density for individual electrodes
        for n_electrodes=1:total_electrodes
            figure('visible', 'off');
            
            new_electrode_number = n_electrodes;
            original_electrode_number = electrodes_order(new_electrode_number); 
    
            plot(spike_density_for_electrode_observation(new_electrode_number, :))
            xlabel('Original Electrode ' + string(original_electrode_number))
            ylabel('Spike Counts binned to 1 second interval')
            filename_png = strcat(file_directory_to_save_data, 'figures/SpikeDensity/', 'Electrode', string(original_electrode_number), '_SpikeDensity.png');
            saveas(gcf, filename_png);
    
            filename_fig = strcat(file_directory_to_save_data, 'figures/SpikeDensity/', 'Electrode', string(original_electrode_number), '_SpikeDensity.fig');
            saveas(gcf, filename_fig);
        end
    end

end

%% THE END