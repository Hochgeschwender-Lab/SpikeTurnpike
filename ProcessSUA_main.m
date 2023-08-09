% Run this after Kilosort and Phy to generate the all_data struct and
% associated metadata.

clear;
w = warning ('off', 'all');

%% Get user inputs
projectFolder = uigetdir('Select the project folder in which you want to analyze multiple groups');
projectFolder = fullfile(projectFolder,'SpikeStuff');

% big UI input window
prompts = {'Cell types set to use for classification (1 for cortical FS/RS, 2 for striatal):',...
          'Timestamps? (0 for none, 1 for whisker stim, 2 for LED stim)',...
          'Stimulation duration in milliseconds (only matters if you have timestamps):',...
          'Cutoff after stim onset to calculate FR (for whole duration, use the value above):'};
default_inputs = {'1', '0', '500', '500'};

user_inputs = inputdlg(prompts, 'User Params', 1, default_inputs);
cell_types_set = str2num(user_inputs{1});
timestamps_type = str2num(user_inputs{2});
stim_duration_ms = str2num(user_inputs{3});
stim_duration_samples = stim_duration_ms*30;
user_time_cutoff_ms = str2num(user_inputs{4});

% cell_types_set = inputdlg('ID of cell types set to use for classification:','Cell Types Set',[1 40],{'1'});
% cell_types_set = str2num(cell_types_set{1});
% 
% timestamps_YorN = inputdlg('Are there any timestamps to extract? (1 for yes, 0 for no)','Timestamps?',[1 40],{'0'});
% timestamps_YorN = str2num(timestamps_YorN{1});
% 
% WhiskerOrOpto = inputdlg('What kind of stimulation did you do? (0 for none, 1 for whisker, 2 for 8 Hz opto)','Stimulation?',[1 40],{'0'});
% WhiskerOrOpto = str2num(WhiskerOrOpto{1});

if timestamps_type == 1
    trialTagsLabels = {'Zero','Low','Mid','Max'};
elseif timestamps_type == 2
    trialTagsLabels = {'8 Hz LED'};
% Add more types of stimulation as elseif conditions here
end

% change stim duration if it's different!
% stim_duration_ms = 500;

% user_time_cutoff_ms = inputdlg('To calculate FRs: how long after the onset of stimulation do you actually want to analyze? (Enter 0 for entire stim duration)','FR time cutoff',[1 40],{'0'});
% user_time_cutoff_ms = str2num(user_time_cutoff_ms{1});

%% Construct Gaussian filter for PSTH smoothing
gausssigma=3;
gausswindow=21;
tempfil = exp(- (-(gausswindow-1)/2:(gausswindow-1)/2)'.^2 ./ (2*gausssigma^2));
tempfilter = tempfil/sum(tempfil); % sum(tempfil) == sqrt(2*pi)*gausssigma, so it's the correct denominator

%% Get electrodes order
%electrodes_order = load(electrode_order_filepath).electrodes_order;
electrodes_order = [14; 20; 16; 18; 1; 31; 3; 29; 5; 27; 7; 25; 9; 23; 11; 21;...
                        13; 19; 15; 17; 12; 22; 10; 24; 8; 26; 6; 28; 4; 30; 2; 32];

%% directory stuff
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
for groupNum = 1:numGroups
    [~,this_group] = fileparts(groupfoldernames(groupNum));
    fprintf('Loading data for %s group\n',this_group);

    groupDir = groupfoldernames{groupNum};
    dinfo2 = dir(groupDir);
    dinfo2(~[dinfo2.isdir]) = []; 
    dinfo2(ismember({dinfo2.name},{'.','..'})) = []; %remove non directories from struct dinfo2
    recfoldernames = fullfile(groupDir, {dinfo2.name}); %create a 1xN cell array for each recording within group, N is number of recordings
    numRecordings = length(recfoldernames);   

    for recNum = 1:numRecordings %loop through the recordings within a group folder
        [~,this_recording] = fileparts(recfoldernames(recNum));
        fprintf('    Loading data for recording %s\n',this_recording);

        recDir = recfoldernames{recNum};

        %% Find SUA, MUA folder with corresponding data 
        SUA_Directory = fullfile(recDir,"SUA"); %path of the MUA directory 
        
        MUA_Directory = fullfile(recDir,"MUA"); %path of the MUA directory  
        
        MUA_allData_Directory = fullfile(recDir,"MUA","allData"); %path of the allData directory
        MUA_figures_Directory = fullfile(recDir,"MUA","figures"); %path of the figures directory

        % Check if figures folder exists within the KS directory folder
        check_if_folder_exists(fullfile(SUA_Directory,'figures','post_sorted'))
        folder_to_save_figure = fullfile(SUA_Directory,'figures','post_sorted');
        
        if timestamps_type ~= 0
            %Extract data needed from MUA folders within a recording of Group N
            timestamp_filepath = fullfile(MUA_allData_Directory,'timestamp_ms.mat');
            % Get timestamps
            try
                timestamp_ms = load(timestamp_filepath).timestamp_ms;
            catch
                timestamp_ms = load(timestamp_filepath);
            end
            % Reject glitched timestamps
            goodTrialInds = (timestamp_ms(:,1) > 0);
            timestamp_ms = timestamp_ms(goodTrialInds,:);
        end

        %electrode_order_filepath = fullfile(MUA_allData_Directory, 'electrodes_order.mat');

        cluster_group_filepath = fullfile(SUA_Directory, 'cluster_group.tsv');
        cluster_info_filepath = fullfile(SUA_Directory, 'cluster_info.tsv');
        
        keepCID_filepath = fullfile(SUA_Directory, 'keepCID.csv');

        % intiate the struct for good_cluster_ids which will be appened for
        % each recording in the project folder
        good_cluster_ids.(this_recording) = [];

        % GRANT: extract good cell IDs from cluster_group.tsv
        cluster_group = tdfread(cluster_group_filepath);
        good_cids = [];
        for ii = 1:length(cluster_group.cluster_id)
            cluster_group_label = strtrim(cluster_group.group(ii,:));
            if strcmp(cluster_group_label, 'good')
                good_cids = [good_cids cluster_group.cluster_id(ii)];
            end
        end

        total_good_ids = length(good_cids);
        
        clear cluster_group

        %% clusterinfo.tsv
        cluster_info = tdfread(cluster_info_filepath);

        all_template_channels = zeros(1, total_good_ids);

        for n_good_channels=1:total_good_ids
            for n_template_channels=1:numel(cluster_info.ch)

                template_channel_phy = cluster_info.ch(n_template_channels);
                template_channel_MATLAB = template_channel_phy+1;

                template_channel_original = find(template_channel_MATLAB == electrodes_order);

                if good_cids(n_good_channels) == cluster_info.cluster_id(n_template_channels)
                    all_template_channels(1, n_good_channels) = template_channel_original;
                end

            end
        end

        clear cluster_info

        %% Load length AND chem stim from files with extension ".ns6" and '.nev' in specified directory of the recording         
        ns6_file_name = dir(fullfile(recDir, '*.ns6'));
        ns6_file_path = fullfile(recDir, ns6_file_name.name);

        % load the correct channel to get an example trace of the
        % stimulator voltage output (only from the first recording)
        if (groupNum == 1)&&(recNum == 1)&&(timestamps_type == 1) % whisker
            ns6_file = openNSx(ns6_file_path, 'read', 'c:35', 'uV', 'precision','double');

            stim_onsets = timestamp_ms(:, 1)'*30; % in samples % 30kHz sampling rate
            %stim_offsets = timestamp_ms(:, 2)'*30; % in samples % 30kHz sampling rate
            stim_intensity = timestamp_ms(:, 3)';

            stim_onsets_low = stim_onsets(stim_intensity == 2);
            stim_onset_low = stim_onsets_low(1);
            StimVoltageTraces_samples(:,1) = ns6_file.Data(stim_onset_low:stim_onset_low+stim_duration_samples)';

            stim_onsets_mid = stim_onsets(stim_intensity == 3);
            stim_onset_mid = stim_onsets_mid(100);
            StimVoltageTraces_samples(:,2) = ns6_file.Data(stim_onset_mid:stim_onset_mid+stim_duration_samples)';
            
            stim_onsets_max = stim_onsets(stim_intensity == 4);
            stim_onset_max = stim_onsets_max(200);
            StimVoltageTraces_samples(:,3) = ns6_file.Data(stim_onset_max:stim_onset_max+stim_duration_samples)';

            save(fullfile(projectFolder,'StimVoltageTraces_samples.mat'), 'StimVoltageTraces_samples');

            % downsample to milliseconds using median of every 30 samples
            tt = array2timetable(StimVoltageTraces_samples, 'SampleRate',30000);
            ttMedian = retime(tt, 'regular', 'median', 'SampleRate',1000);
            StimVoltageTraces_ms = [ttMedian.StimVoltageTraces_samples1 ttMedian.StimVoltageTraces_samples2 ttMedian.StimVoltageTraces_samples3];
            save(fullfile(projectFolder,'StimVoltageTraces_ms.mat'), 'StimVoltageTraces_ms');
            clear tt ttMedian
        elseif (groupNum == 1)&&(recNum == 1)&&(timestamps_type == 2) % LED/opto
            ns6_file = openNSx(ns6_file_path, 'read', 'c:33', 'uV', 'precision','double');

            stim_onsets = timestamp_ms(:, 1)'*30; % in samples % 30kHz sampling rate

            stim_onset = stim_onsets(10);
            % startSample = stim_onset - stim_duration_samples;
            % endSample = stim_onset + 2*stim_duration_samples - 1;

            %StimVoltageTrace_samples = ns6_file.Data(startSample:endSample)';
            StimVoltageTrace_samples = ns6_file.Data(stim_onset:stim_onset+stim_duration_samples)';
            save(fullfile(projectFolder,'StimVoltageTrace_samples.mat'), 'StimVoltageTrace_samples');

            % downsample to milliseconds using median of every 30 samples
            tt = timetable(StimVoltageTrace_samples, 'SampleRate',30000);
            ttMedian = retime(tt, 'regular', 'median', 'SampleRate',1000);
            StimVoltageTrace_ms = ttMedian.StimVoltageTrace_samples;
            save(fullfile(projectFolder,'StimVoltageTrace_ms.mat'), 'StimVoltageTrace_ms');
            clear tt ttMedian
        else
            ns6_file = openNSx(ns6_file_path, 'c:1'); % doesn't matter for chem stim
        end

        duration = ns6_file.MetaTags.DataDurationSec; %load length of recording

        nev_file_name = dir(fullfile(recDir, '*.nev'));
        nev_file_path = fullfile(recDir, nev_file_name.name);


        % make this into an option 
           
        % %pull out when a chemical stimulation occurred, in this case,
        % %assumes one comment for when a treatment occurred 
        % nev_file = openNEV(nev_file_path);
        % chem_stim_time_s = nev_file.Data.Comments.TimeStampStartedSec;
        % chem_stim_time_samples = nev_file.Data.Comments.TimeStampStarted;
        % chem_stim_time_note = nev_file.Data.Comments.Text;


        clear ns6_file nev_file

        %% Create struct to organize data
        for n_fieldnames=1:total_good_ids
            cid = good_cids(n_fieldnames);

            fieldname = strcat('cid', string(cid));

            
            good_cluster_ids.(this_recording).(fieldname) = struct('Header', [], 'Template_Channel', [], 'Template_Channel_Position', [], "SpikeTimes_all", [], "SpikeTimes_baseline", [], "SpikeTimes_stim", [],...
                'Recording_Duration', duration, 'Cell_Type', [], 'SpikeTrains_trials', [], 'SpikeTrains_baseline', [], 'SpikeTrains_stim', [], 'FRs_baseline', [], 'FRs_stim', [], 'MeanFR_baseline', [], 'MeanFR_stim', [],...
                'MeanFR_total', [], 'Mean_Waveform', [],'UnNormalized_Template_Waveform', [], 'Normalized_Template_Waveform', [], 'TroughToPeak_duration', [], 'Peak2ToTrough_ratio', [], 'PeakToPeak_ratio', [],...
                'Peak1ToTrough_ratio', [], 'SpikeHalfWidth', [], 'Amplitude', [], 'Stim_Onsets_samples', [], 'peak1_normalized_amplitude', [], 'ISI_violations_percent', [], 'Stim_Offsets_samples', [],...
                'Stim_Intensity', [], 'Sampling_Frequency', 30000, 'SpikeTimes_trials', [], 'StimProb', [], 'StimResponsivity', [], 'IsSingleUnit', [], 'SpikeTrains_for_PSTHs', [], 'PSTHs_raw', [], 'PSTHs_conv', [],...
                'FR_time_cutoff_after_stim_ms', user_time_cutoff_ms, 'PeakEvokedFR', [], 'PeakEvokedFR_Latency', [], 'ModulationIndex', [], 'MeanFR_inst_baseline', [], 'MeanFR_inst_stim', [], ...
                'FirstSpikeLatency_perTrial', [], 'FRs_baseline_vec', [], 'ISI_baseline_vec', [], 'ISI_baseline_CV', [], 'FanoFactor_baseline', [], 'FanoFactor_stim', [],...
                'ISI_pdf_x',[], 'ISI_pdf_y',[], 'ISI_pdf_peak_xy',[], 'FirstSpikeLatency_pdf_y',[], 'FirstSpikeLatency_pdf_x',[], 'FirstSpikeLatency', [], 'FirstSpikeLatency_Reliability',[],...
                'SpikeTrains_baseline_ms',[], 'SpikeTrains_trials_ms', [], 'SpikeTrains_stim_ms', [], 'DeflectionResponseRatios',[], 'Depth',[]);
            % NOTE: relative_time_ms and relative_time_samples are the time
            % arrays for the spike trains and PSTHs, stored in SpikeStuff
            % with all the other global metadata.
            
            % make this into an option
            % good_cluster_ids.(this_recording).(fieldname) = struct('Header', [], 'Template_Channel', [], 'Template_Channel_Position', [], "SpikeTimes_all", [], "SpikeTimes_baseline", [], "SpikeTimes_stim", [],...
            %     'Recording_Duration', duration, 'Cell_Type', [], 'SpikeTrains_trials', [], 'SpikeTrains_baseline', [], 'SpikeTrains_stim', [], 'FRs_baseline', [], 'FRs_stim', [], 'MeanFR_baseline', [], 'MeanFR_stim', [],...
            %     'MeanFR_total', [], 'Mean_Waveform', [],'UnNormalized_Template_Waveform', [], 'Normalized_Template_Waveform', [], 'TroughToPeak_duration', [], 'Peak2ToTrough_ratio', [], 'PeakToPeak_ratio', [],...
            %     'Peak1ToTrough_ratio', [], 'SpikeHalfWidth', [], 'Amplitude', [], 'Stim_Onsets_samples', [], 'peak1_normalized_amplitude', [], 'ISI_violations_percent', [], 'Stim_Offsets_samples', [],...
            %     'Stim_Intensity', [], 'Sampling_Frequency', 30000, 'SpikeTimes_trials', [], 'StimProb', [], 'StimResponsivity', [], 'IsSingleUnit', [], 'ChemStimTime_s', chem_stim_time_s, 'ChemStimTime_samples', chem_stim_time_samples,...
            %         'ChemStimTime_note',chem_stim_time_note);

            % Arrange the fieldnames in alphabetical order
            [~, alphabetical_order] = sort(fieldnames(good_cluster_ids.(this_recording).(fieldname)));
            good_cluster_ids.(this_recording).(fieldname) = orderfields(good_cluster_ids.(this_recording).(fieldname));

            % This field is only here to add special notes
            %header = {'all fields are in samples except trough to peak which is in ms'};
            %good_cluster_ids.(this_recording).(fieldname).Header = header;

        end

        %% Load channel_positions (xy coordinates)
        channel_positions_filepath = fullfile(SUA_Directory, 'channel_positions.npy');
        channel_positions = readNPY(channel_positions_filepath);

        %% Run getWaveFroms.m

        sp = loadKSdir(SUA_Directory);

        %FRAC = 0.008;%0.02%0.015;%0.0002; % fraction of all spikes to sample from
        gwfparams.sr = sp.sample_rate;
        gwfparams.dataDir = SUA_Directory;
        gwfparams.fileName = 'temp_wh.dat';
        gwfparams.dataType = sp.dtype;
        gwfparams.nCh = sp.n_channels_dat;
        gwfparams.wfWin = [-(0.002*gwfparams.sr) 0.002*gwfparams.sr];              % Number of samples before and after spiketime to include in waveform
        %gwfparams.nWf = floor(length(sp.st) * FRAC);
        if length(sp.st) >= 2000 % use 2000 waveforms at most
            gwfparams.nWf = 2000;
        else
            gwfparams.nWf = length(sp.st);
        end
        gwfparams.spikeTimes = ceil(sp.st * sp.sample_rate);
        gwfparams.spikeClusters = sp.clu;

        clear sp

        wf = getWaveForms(gwfparams);

        all_spike_clusters = gwfparams.spikeClusters;
        all_spike_times = gwfparams.spikeTimes; % in samples
        all_unitIDs = wf.unitIDs;
        all_waveforms = wf.waveForms;
        all_meanwaveforms = wf.waveFormsMean;

        [total_unitIDs, total_nwf, total_electrodes, total_samples] = size(all_waveforms);

        clear gwfparams wf

        %% Get the spike_cluster and the spike_time from gwfparams

        for n_clustertimes=1:length(all_spike_clusters)

            is_this_a_good_cid = find(all_spike_clusters(n_clustertimes) == good_cids);

            if isempty(is_this_a_good_cid)
                continue;

            else

                cid = good_cids(is_this_a_good_cid);
                fieldname = strcat('cid', string(cid));

                spike_time = all_spike_times(n_clustertimes);

                good_cluster_ids.(this_recording).(fieldname).SpikeTimes_all =...
                    [good_cluster_ids.(this_recording).(fieldname).SpikeTimes_all spike_time]; % in samples % 30kHz sampling rate

            end

        end

        clear all_spike_clusters all_spike_times

        %% After filtering - Not normalized
        %% Bandpass filtered data
        low_frequency_pass = 300; % Hz % From Manny's code.
        low_frequency_stop = 200; % Hz % From Prakash et al.

        % Define high cutoff frequency
        high_frequency_stop = 7000; % Hz % From Manny's code.
        high_frequency_pass = 6000; % Hz % From Prakash et al.

        % Define sampling_frequency
        sampling_frequency = 30e3; % Hz % Downsampled from 30kHz

        % Find nyquist_frequency for normalizing the frequency
        nyquist_frequency = (sampling_frequency/2);

        % Determine passband and stopband edge frequency
        Wp = [low_frequency_pass high_frequency_pass] / nyquist_frequency;
        Ws = [low_frequency_stop high_frequency_stop] / nyquist_frequency;

        % Assign dB of ripple and attenuation
        [N, Wn] = buttord(Wp, Ws, 3, 20);
        [B, A] = butter(N, Wn);

        filtered_waveform = zeros(total_good_ids, total_samples);

        figure('visible', 'off');
        % figure;
        tiledlayout('flow');

        for n_cell=1:total_good_ids

            nexttile;

            selected_cid = good_cids(n_cell);
            fieldname = strcat('cid', string(selected_cid));

            index_in_unitID = find(selected_cid == all_unitIDs);
            %interested_template = all_template_channels(n_cell);

            particular_waveform = squeeze(all_waveforms(index_in_unitID, :, :, :)); % [waveforms x channels x samples] for this cluster
            mean_across_nwf = squeeze(nanmean(particular_waveform, 1)); % mean across individual waveforms: now [channels x samples]

            good_cluster_ids.(this_recording).(fieldname).Mean_Waveform =...
                [good_cluster_ids.(this_recording).(fieldname).Mean_Waveform mean_across_nwf]; % in samples % 30kHz sampling rate

            [interested_template,~] = find(mean_across_nwf==min(min(mean_across_nwf))); % GRANT: find channel where the mean template waveform has the greatest amplitude
            
            unit_amplitude = abs(min(min(mean_across_nwf))); % getting and saving mean amplitude (absolute value)
            good_cluster_ids.(this_recording).(fieldname).Amplitude =...
                [good_cluster_ids.(this_recording).(fieldname).Amplitude unit_amplitude];

            % GRANT: this gives the actual channel that the waveform occurs on (Phy is sometimes wrong)
            good_cluster_ids.(this_recording).(fieldname).Template_Channel = ...
                [good_cluster_ids.(this_recording).(fieldname).Template_Channel interested_template];

            % Also save the channel position (X/Y) in micrometers
            channel_position = channel_positions(interested_template,:);
            good_cluster_ids.(this_recording).(fieldname).Template_Channel_Position = ...
                [good_cluster_ids.(this_recording).(fieldname).Template_Channel_Position channel_position];
            good_cluster_ids.(this_recording).(fieldname).Depth = channel_position(2);

            %interested_waveform = mean_across_nwf(interested_template, :);
            interested_waveform = squeeze(all_meanwaveforms(index_in_unitID,interested_template,:));

            good_cluster_ids.(this_recording).(fieldname).UnNormalized_Template_Waveform =...
                [good_cluster_ids.(this_recording).(fieldname).UnNormalized_Template_Waveform interested_waveform]; % in samples % 30kHz sampling rate

            filtered_waveform(n_cell, :) = filtfilt(B, A, interested_waveform);
            %filtered_waveform(n_cell, :) = interested_waveform;

            plot(filtered_waveform(n_cell, :))
            xlabel('cid' + string(selected_cid))
            hold on;

        end
        hold off;
        set(gcf, 'Position', get(0, 'Screensize'));
        saveas(gcf, fullfile(folder_to_save_figure, 'all_cell_af_fil_unnormalized.png'))

        clear all_unitIDs all_waveforms all_meanwaveforms

        %% After filtering - Normalized
        figure('visible', 'off');
        tiledlayout('flow')

        for n_cell=1:total_good_ids
            nexttile;

            selected_cid = good_cids(n_cell);
            fieldname = strcat('cid', string(selected_cid));

            interested_waveform = filtered_waveform(n_cell, :);
            normalized_waveform = interested_waveform/abs(min(interested_waveform));

            good_cluster_ids.(this_recording).(fieldname).Normalized_Template_Waveform =...
                [good_cluster_ids.(this_recording).(fieldname).Normalized_Template_Waveform normalized_waveform']; % in samples % 30kHz sampling rate

            % Normalizing
            trough_location = find(min(normalized_waveform) == normalized_waveform);
            max_after_trough = max(normalized_waveform(trough_location:end));
            max_after_trough_location = find(max_after_trough == normalized_waveform);
            max_before_trough = max(normalized_waveform(1:trough_location));
            max_before_trough_location = find(max_before_trough == normalized_waveform);
            difference = (max_after_trough_location - trough_location)/30; % ms trough to peak
            %difference = (max_after_trough_location - max_before_trough_location)/30; % ms peak to peak

            % save trough to peak duration (milliseconds)
            good_cluster_ids.(this_recording).(fieldname).TroughToPeak_duration = ...
                [good_cluster_ids.(this_recording).(fieldname).TroughToPeak_duration difference];

            % peak2-trough ratio
            trough_value = normalized_waveform(trough_location);
            trough_abs_value = abs(trough_value);
            peak2_value = abs(normalized_waveform(max_after_trough_location));
            peak2ToTroughRatio = peak2_value / trough_abs_value;
            good_cluster_ids.(this_recording).(fieldname).Peak2ToTrough_ratio = ...
                [good_cluster_ids.(this_recording).(fieldname).Peak2ToTrough_ratio peak2ToTroughRatio];

            % peak1 amplitude and peak-peak ratio
            peak1_value = normalized_waveform(max_before_trough_location);
            good_cluster_ids.(this_recording).(fieldname).peak1_normalized_amplitude = ...
                [good_cluster_ids.(this_recording).(fieldname).peak1_normalized_amplitude peak1_value];

            peak1_abs_value = abs(peak1_value);
            peak2_value = abs(normalized_waveform(max_after_trough_location));
            peakToPeak_ratio = peak1_abs_value / peak2_value;
            good_cluster_ids.(this_recording).(fieldname).PeakToPeak_ratio = ...
                [good_cluster_ids.(this_recording).(fieldname).PeakToPeak_ratio peakToPeak_ratio];

            % peak1-trough ratio
            peak1ToTroughRatio = peak1_abs_value / trough_abs_value;
            good_cluster_ids.(this_recording).(fieldname).Peak1ToTrough_ratio = ...
                [good_cluster_ids.(this_recording).(fieldname).Peak1ToTrough_ratio peak1ToTroughRatio];

            % spike half-width
            halfAmplitude = normalized_waveform(trough_location) / 2;

            wf_peak1_to_trough = normalized_waveform(max_before_trough_location:trough_location); % only looking at the depolarization slope
            [~,ind] = min(abs(wf_peak1_to_trough - halfAmplitude));
            value1 = wf_peak1_to_trough(ind);
            inds1 = find(normalized_waveform == value1);
            ind1 = inds1(end);

            wf_trough_to_peak2 = normalized_waveform(trough_location:max_after_trough_location); % only looking at the repolarization slope
            [~,ind] = min(abs(wf_trough_to_peak2 - halfAmplitude));
            value2 = wf_trough_to_peak2(ind);
            inds2 = find(normalized_waveform == value2);
            ind2 = inds2(1);

            halfWidth = (ind2-ind1)/30;
            good_cluster_ids.(this_recording).(fieldname).SpikeHalfWidth = ...
                [good_cluster_ids.(this_recording).(fieldname).SpikeHalfWidth halfWidth];

            % Quality control: proportion of ISI violations (counted as any < 2 ms) should be < 0.01%
            % Calculate PROP(isi > 2s), the proportion of inter-spike
            % intervals that exceed two seconds
            spike_times = good_cluster_ids.(this_recording).(fieldname).SpikeTimes_all;
            isi_vec = zeros(length(spike_times)-1,1);
            for spike_ind = 1:length(spike_times)-1
                isi_vec(spike_ind) = spike_times(spike_ind+1) - spike_times(spike_ind);
            end
            isi_vec_ms = isi_vec/30; % conversion from samples (assuming 30 kHz sampling rate) to milliseconds

            ISI_violations_count = sum(isi_vec_ms < 2);
            ISI_violations_percent = (ISI_violations_count / length(isi_vec_ms)) * 100;
            good_cluster_ids.(this_recording).(fieldname).ISI_violations_percent = ISI_violations_percent;

            % if continuous data (no trials), calculate ISI stuff here and
            % store ISI vector
            if timestamps_type == 0
                good_cluster_ids.(this_recording).(fieldname).ISI_baseline_vec = isi_vec_ms;
                % calculate and store baseline ISI probability density function (PDF)
                [ISI_pdf_x, ISI_pdf_y] = EmpiricalPDF(isi_vec_ms);
                good_cluster_ids.(this_recording).(fieldname).ISI_pdf_y = ISI_pdf_y;
                good_cluster_ids.(this_recording).(fieldname).ISI_pdf_x = ISI_pdf_x;
                % ... peak of the baseline ISI PDF ("ISI Peak" in Shin & Moore, 2019)
                [ISI_pdf_peak_y, ISI_pdf_peak_ind] = max(ISI_pdf_y);
                ISI_pdf_peak_x = ISI_pdf_x(ISI_pdf_peak_ind);
                good_cluster_ids.(this_recording).(fieldname).ISI_pdf_peak_xy = [ISI_pdf_peak_x, ISI_pdf_peak_y];
            
                % ... coefficient of variation (CV) of ISIs
                good_cluster_ids.(this_recording).(fieldname).ISI_baseline_CV = std(isi_vec_ms) / mean(isi_vec_ms);
            end

            if cell_types_set == 1
                %% Classify cortical RS and FS
%                 cell_types = {'RS','FS','Unclassified'};
%                 if (halfWidth <= 0.2) && (ISI_violations_percent < 0.01)
%                     good_cluster_ids.(file_number_name).(fieldname).Cell_Type = 'FS';
%                 elseif (halfWidth > 0.2) && (ISI_violations_percent < 0.01)
%                     good_cluster_ids.(file_number_name).(fieldname).Cell_Type = 'RS';
%                 else
%                     good_cluster_ids.(file_number_name).(fieldname).Cell_Type = 'Unclassified';
%                 end

                cell_types = {'RS','FS'};
                if difference <= 0.4  %was 0.5 before 5.12.23 changed
                    good_cluster_ids.(this_recording).(fieldname).Cell_Type = 'FS';
                else
                    good_cluster_ids.(this_recording).(fieldname).Cell_Type = 'RS';
                end
            end

            if cell_types_set == 2
                %% GRANT: classifications for HD striatum cell types from Yamin
                %        et al. (2021), excluding the "baseline FR" parameter
    
                % Valley to peak duration (VPD) is already given by the
                % variable "difference"
                
                % Calculate overall FR and save as Total_Avg_FR
                overall_FR = length(spike_times)/duration; % calculate FR in spikes per second
                good_cluster_ids.(this_recording).(fieldname).MeanFR_total = overall_FR;
           
                % Calculate PROP(isi > 2s), the proportion of inter-spike
                % intervals that exceed two seconds
                isi_vec_s = isi_vec/30000; % conversion from samples (assuming 30 kHz sampling rate) to seconds
                num_longISIs = sum(isi_vec_s > 2);
                prop_LongISIs = num_longISIs / length(isi_vec_s);
                good_cluster_ids.(this_recording).(fieldname).prop_LongISIs = prop_LongISIs;
    
                % Classify the cell
                cell_types = {'MSN','FSI','UIN','TAN','Other'};
                if (difference > 0.5) && (prop_LongISIs >= 0.4)
                    good_cluster_ids.(this_recording).(fieldname).Cell_Type = 'MSN'; % medium spiny neuron
                elseif (difference <= 0.5) && (prop_LongISIs < 0.3)
                    good_cluster_ids.(this_recording).(fieldname).Cell_Type = 'FSI'; % fast-spiking interneuron
                elseif (difference <= 0.5) && (prop_LongISIs > 0.4) && (overall_FR < 4)
                    good_cluster_ids.(this_recording).(fieldname).Cell_Type = 'UIN'; % unidentified neuron (spooky! Area 51!)
                elseif (difference > 0.5) && (prop_LongISIs < 0.4)
                    good_cluster_ids.(this_recording).(fieldname).Cell_Type = 'TAN'; % tonically active neuron
                else
                    good_cluster_ids.(this_recording).(fieldname).Cell_Type = 'Other'; % just so every cell has a "type"
                end

                % cell_types = {'MSN','FSI'};
                % if difference < 0.35
                %     good_cluster_ids.(this_recording).(fieldname).Cell_Type = 'FSI';
                % else
                %     good_cluster_ids.(this_recording).(fieldname).Cell_Type = 'MSN';
                % end

            end

            %% Quality control: classify as MUA or SUA based on ISI violations
            if ISI_violations_percent < 0.1
                good_cluster_ids.(this_recording).(fieldname).IsSingleUnit = 1;
            else
                good_cluster_ids.(this_recording).(fieldname).IsSingleUnit = 0;
            end

            plot(normalized_waveform)
            xlabel('cid' + string(selected_cid))
            hold on;
%             xline(trough_location, '--r', {'Trough'})
%             xline(max_after_trough_location, '--r', {'Peak after Trough'})
            title('Difference ' + string(round(difference, 2)) + ' ms, ' + good_cluster_ids.(this_recording).(fieldname).Cell_Type)
            plot(max_before_trough_location,max_before_trough, '.', 'Color','r');
            plot(trough_location,trough_value, '.', 'Color','r');
            plot(max_after_trough_location,max_after_trough, '.', 'Color','r');
            line([ind1 ind2],[value1 value2], 'Color','k');


        end
        hold off;
        set(gcf, 'Position', get(0, 'Screensize'));
        saveas(gcf, fullfile(folder_to_save_figure, 'all_cell_af_fil_normalized.png'))

        clear filtered_waveform normalized_waveform interested_waveform

        %% Do FR stuff
        if timestamps_type ~= 0 % GRANT: want to skip this whole section if there are no timestamps for manipulations

            % Get single-recording dataset
            rec_data = good_cluster_ids.(this_recording);

            % Run Grant's cool new function to do all the FR stuff, then
            % plug rec_data back into the main dataset
            relative_time_ms = -stim_duration_ms:2*stim_duration_ms-1; % because indexing (in CalculateEvokedFRs_OneRecording) is inclusive of the trial start and exclusive of the trial end
            % ^ Time 0 is on index 501 because MATLAB indexes from 1 and
            % all of the spike times are translated from timestamps to
            % indices in the same way (in CalculateEvokedFRs_OneRecording)
            % to avoid the 'indexing from 0' error.
            % SO, considering this, every spike or timestamp is at INDEX
            % t+1 if TIME t is where it actually falls.
            
            relative_time_samples = -stim_duration_samples:2*stim_duration_samples-1;

            rec_data = CalculateEvokedFRs_OneRecording(rec_data,timestamp_ms,stim_duration_ms,user_time_cutoff_ms,tempfilter,relative_time_ms);
            
            good_cluster_ids.(this_recording) = rec_data;
        end
        
        %% TO DO: if stim was optogenetic, find opto-tagged units
        if timestamps_type == 2

        end

    end

    %% put single-group data into all_data
    all_data.(this_group) = good_cluster_ids;

    clear good_cluster_ids

end

%% Save the struct data
save(fullfile(projectFolder, 'all_data.mat'), 'all_data', '-v7.3');
save(fullfile(projectFolder, 'channel_positions.mat'), 'channel_positions');
save(fullfile(projectFolder, 'electrodes_order.mat'), 'electrodes_order');
save(fullfile(projectFolder, 'cell_types.mat'), 'cell_types');
if timestamps_type ~= 0
    save(fullfile(projectFolder, 'trialTagsLabels.mat'), 'trialTagsLabels');
    save(fullfile(projectFolder, 'relative_time_ms.mat'), 'relative_time_ms');
    save(fullfile(projectFolder, 'relative_time_samples.mat'), 'relative_time_samples');
end