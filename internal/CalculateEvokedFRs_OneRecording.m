function rec_data = CalculateEvokedFRs_OneRecording(rec_data,timestamp_ms,stim_duration_ms,user_time_cutoff_ms,tempfilter,relative_time_ms)
% That big part of The Shashwat, simplified. Calculates baseline and evoked
% firing rates per trial for a single unit/cell. Also stores spike trains
% for each. Among many other things.
%
% INPUTS
% - rec_data: one recording from the struct all_data, e.g.,
%       all_data.Control.(Recording_Name)
% - timestamp_ms: [numTrials x 3] array where column 1 is stimulation
%       onsets, column 2 is stimulation offsets, and column 3 is 
%       stimulation intensity/type. Columns 1-2 are in milliseconds. Should
%       correspond to the specific recording.
% - stim_duration_ms: duration of the stimulation window in ms. Also
%       determines the duration of the baseline window and post-stimulation
%       window.
% - user_time_cutoff_ms: milliseconds after the stimulation onset to calculate
%       FRs in. Set as 0 to measure FR for the entire stim. duration.
%       E.g., if you want to get FRs during the first 100 ms of
%       stimulation, set user_time_cutoff to 100.
% - tempfilter: Gaussian filter defined elsewhere. For example,
%       gausssigma=3;
%       gausswindow=21;
%       tempfil = exp(- (-(gausswindow-1)/2:(gausswindow-1)/2)'.^2 ./ (2*gausssigma^2));
%       tempfilter = tempfil/sum(tempfil); % sum(tempfil) == sqrt(2*pi)*gausssigma, so it's the correct denominator
% - chem_timestamp: sample number where a chemical (e.g., CTZ) was added.

stim_duration_samples = stim_duration_ms*30;

stim_onsets = timestamp_ms(:, 1)'*30; % in samples % 30kHz sampling rate
stim_offsets = timestamp_ms(:, 2)'*30; % in samples % 30kHz sampling rate
stim_intensity = timestamp_ms(:, 3)';

numTrials = length(stim_onsets);

% accounting for varying stim. intensities/types
uniqueTrialTags = unique(stim_intensity);

cellIDs = fieldnames(rec_data);
for cellID_num = 1:length(cellIDs) % loop through units in recording
    cellID = cellIDs{cellID_num};
    
    spikeTimes_all = rec_data.(cellID).SpikeTimes_all; % array of spike timestamps (in samples) from entire recording

    % storing onset/offset/intensity metadata per unit
    rec_data.(cellID).Stim_Onsets_samples = stim_onsets; % formerly Stimulation_Onsets
    rec_data.(cellID).Stim_Offsets_samples = stim_offsets; % formerly Stimulation_Offsets
    rec_data.(cellID).Stim_Intensity = stim_intensity; % formerly Stimulation_Intensity

    spikeTimes_trials = cell(1,numTrials); % this would be "evoked_spike_times_particular_trial" in The Shashwat
    spikeTimes_baseline = cell(1,numTrials); % "baseline_spike_times_particular_trial"
    spikeTimes_stim = cell(1,numTrials); % "during_stim_spike_times_particular_trial"

    % spikeTrains_trials = cell(1,numTrials); % this would be "evoked_spike_particular_trial" in The Shashwat
    % spikeTrains_baseline = cell(1,numTrials); % "baseline_spike_particular_trial"
    % spikeTrains_stim = cell(1,numTrials);
    spikeTrains_trials = int8(zeros(numTrials,3*stim_duration_samples));
    spikeTrains_baseline = int8(zeros(numTrials,stim_duration_samples)); % changed these to be [n_trials x n_samples] arrays (5/25/2023) 
    spikeTrains_stim = int8(zeros(numTrials,stim_duration_samples));

    spikeTrains_trials_ms = int8(zeros(numTrials,3*stim_duration_ms));
    spikeTrains_baseline_ms = int8(zeros(numTrials,stim_duration_ms)); % also storing in millisecond inds
    spikeTrains_stim_ms = int8(zeros(numTrials,stim_duration_ms));

    FRs_baseline_vec = zeros(numTrials,1);
    FRs_baseline = cell(1,length(uniqueTrialTags));
    % FRs_stim = zeros(1,numTrials);
    FRs_stim = cell(1,length(uniqueTrialTags));

    % raw PSTHs per trial type (i.e., stim intensity)
    SpikeTrains_for_PSTHs = cell(1,length(uniqueTrialTags));
    FirstSpikeLatency_perTrial = cell(1,length(uniqueTrialTags));

    ISI_baseline_vec = [];

    for trialNum = 1:numTrials

        trialTag = stim_intensity(trialNum);
        trialTagInd = find(uniqueTrialTags==trialTag);

        onset = stim_onsets(trialNum);
        offset = onset + stim_duration_samples;
        trialStart = onset - stim_duration_samples;
        trialEnd = offset + stim_duration_samples;

        % Check if trial window is the proper length
        if trialEnd-trialStart ~= 3*stim_duration_samples
            error("trial length is wack!");
        end

        spikeTrain_trial = int8(zeros(1,3*stim_duration_samples)); % 0 for no spike, 1 for spike. Indexed by sample.
        % note that stim_duration_samples is also the 'relative onset' when
        % indexing into spikeTrain_trial, and 2*stim_duration_samples is
        % the 'relative offset'

        spikeTimes_trial = spikeTimes_all((spikeTimes_all >= trialStart)&(spikeTimes_all < trialEnd)); % 5/15/23: changed from (spikeTimes_all <= trialEnd) because if a spike landed on trialEnd, spikeTrain_trial would be 1 index too long

        % store spike timestamp vectors in cell arrays
        spikeTimes_trials{1,trialNum} = spikeTimes_trial;

        spikeTimes_baseline_trial = spikeTimes_all((spikeTimes_all >= trialStart)&(spikeTimes_all < onset));
        spikeTimes_baseline{1,trialNum} = spikeTimes_baseline_trial;

        spikeTimes_stim_trial = spikeTimes_all((spikeTimes_all >= onset)&(spikeTimes_all < offset));
        spikeTimes_stim{1,trialNum} = spikeTimes_stim_trial;

        % get baseline/spontaneous ISIs if there are at least two spikes
        if length(spikeTimes_baseline_trial) > 1
            ISI_baseline_trial_vec = zeros(length(spikeTimes_baseline_trial)-1,1);
            for ii = 1:length(spikeTimes_baseline_trial)-1
                ISI_baseline_trial_vec(ii,1) = spikeTimes_baseline_trial(ii+1) - spikeTimes_baseline_trial(ii);
            end
            ISI_baseline_trial_vec = ISI_baseline_trial_vec/30; % conversion from samples (assuming 30 kHz sampling rate) to milliseconds
            ISI_baseline_vec = [ISI_baseline_vec; ISI_baseline_trial_vec];
        end

        % convert spikeTimes to spikeTrains, 0s and 1s indexed by time
        spikeInds = spikeTimes_trial - trialStart + 1; % +1 to account for MATLAB indexing from 1. Hence, t = index - 1.

        spikeTrain_trial(spikeInds) = 1; % the Spike Train is boarding at the station

        spikeTrain_baseline = spikeTrain_trial(1:stim_duration_samples);
        spikeTrain_stim = spikeTrain_trial(stim_duration_samples+1:2*stim_duration_samples);

        % store spike trains in cell arrays
        % spikeTrains_trials{1,trialNum} = spikeTrain_trial;
        % spikeTrains_baseline{1,trialNum} = spikeTrain_baseline;
        % spikeTrains_stim{1,trialNum} = spikeTrain_stim;
        spikeTrains_trials(trialNum,:) = spikeTrain_trial;
        spikeTrains_baseline(trialNum,:) = spikeTrain_baseline;
        spikeTrains_stim(trialNum,:) = spikeTrain_stim;

        % store spike trains in the raw PSTH format - basically just
        % organized into a matrix per trial type/tag/label
        spikeTrain_trial_ms = int8(zeros(1,length(spikeTrain_trial)/30));
        for ms_ind = 1:length(spikeTrain_trial_ms)
            spikeTrain_one_ms = spikeTrain_trial(((ms_ind-1)*30)+1:ms_ind*30);
            if sum(spikeTrain_one_ms) >= 1
                spikeTrain_trial_ms(1,ms_ind) = 1;
            end
        end
        SpikeTrains_for_PSTHs{1,trialTagInd} = [SpikeTrains_for_PSTHs{1,trialTagInd}; spikeTrain_trial_ms];
        
        spikeTrains_trials_ms(trialNum,:) = spikeTrain_trial_ms;
        spikeTrains_baseline_ms(trialNum,:) = spikeTrain_trial_ms(1:stim_duration_ms);
        spikeTrains_stim_ms(trialNum,:) = spikeTrain_trial_ms(stim_duration_ms+1:2*stim_duration_ms);

        % get and store 1st spike latency for every trial in the same
        % format as the raw PSTHs are stored above
        first_spike_ind = find(spikeTrain_stim, 1, "first");
        if isempty(first_spike_ind) % no spikes occurred
            FirstSpikeLatency_perTrial{1,trialTagInd} = [FirstSpikeLatency_perTrial{1,trialTagInd}; NaN];
        else % if there was at least one spike
            first_spike_latency_samples = first_spike_ind - 1; % conversion from index to sample time. Index 1 is the onset (inclusive), = time 0.
            first_spike_latency_ms = first_spike_latency_samples/30;
            FirstSpikeLatency_perTrial{1,trialTagInd} = [FirstSpikeLatency_perTrial{1,trialTagInd}; first_spike_latency_ms];
        end

        % calculate and store firing rates
        spikeCount_baseline = sum(spikeTrain_baseline);
        FR_baseline = spikeCount_baseline*(30000/length(spikeTrain_baseline));
        FRs_baseline_vec(trialNum,1) = FR_baseline;
        FRs_baseline{1,trialTagInd} = [FRs_baseline{1,trialTagInd}; FR_baseline];

        if user_time_cutoff_ms == stim_duration_ms % count spikes in entire stim window
            spikeCount_stim = sum(spikeTrain_stim);
            FR_stim = spikeCount_stim*(30000/length(spikeTrain_stim));
        else % count spikes from stim onset to user-defined cutoff
            spikeTrain_stim_userWindow = spikeTrain_stim(1:(user_time_cutoff_ms*30)+1); % being inclusive of the user-defined cutoff
            spikeCount_stim = sum(spikeTrain_stim_userWindow);
            FR_stim = spikeCount_stim*(30000/length(spikeTrain_stim_userWindow));
        end
        %FRs_stim(trialNum) = FR_stim;
        FRs_stim{1,trialTagInd} = [FRs_stim{1,trialTagInd}; FR_stim];

    end

    % store trial-by-trial spike timestamp vectors
    rec_data.(cellID).SpikeTimes_trials = spikeTimes_trials; % formerly "Spike_Times_Evoked"
    rec_data.(cellID).SpikeTimes_baseline = spikeTimes_baseline; % formerly "Spike_Times_Baseline"
    rec_data.(cellID).SpikeTimes_stim = spikeTimes_stim; % formerly "Spike_Times_During"

    % store trial-by-trial spike trains
    rec_data.(cellID).SpikeTrains_trials = spikeTrains_trials; % formerly "Evoked_Spike_Trials"
    rec_data.(cellID).SpikeTrains_baseline = spikeTrains_baseline; % formerly "Baseline_Spike_Trials"
    rec_data.(cellID).SpikeTrains_stim = spikeTrains_stim;
    rec_data.(cellID).SpikeTrains_trials_ms = spikeTrains_trials_ms;
    rec_data.(cellID).SpikeTrains_baseline_ms = spikeTrains_baseline_ms;
    rec_data.(cellID).SpikeTrains_stim_ms = spikeTrains_stim_ms;

    % store trial-by-trial FR vectors
    rec_data.(cellID).FRs_baseline_vec = FRs_baseline_vec;
    rec_data.(cellID).FRs_baseline = FRs_baseline;
    rec_data.(cellID).FRs_stim = FRs_stim; % this field was called "FR_Evoked_Across_Trials" in The Shashwat

    % store trial-by-trial 1st spike latencies
    rec_data.(cellID).FirstSpikeLatency_perTrial = FirstSpikeLatency_perTrial;

    % calculate and store mean baseline FR from all trials (irrespective of
    % tag/type)
    MeanFR_baseline = mean(FRs_baseline_vec);
    rec_data.(cellID).MeanFR_baseline = MeanFR_baseline; % formerly Baseline_Avg_FR
    
    % for trialTagInd = 1:length(uniqueTrialTags)
    %     trialTag = uniqueTrialTags(trialTagInd);
    %     meanFR_stim(trialTagInd,1) = mean(FRs_stim(find(stim_intensity==trialTag)));
    % end
    
    % rec_data.(cellID).MeanFR_stim = MeanFR_stim; % vector indexed by stim types/intensities, which are saved in trialTagsLabels.mat.
    %     % formerly Evoked_Avg_FR

    % store spike trains for PSTHs
    rec_data.(cellID).SpikeTrains_for_PSTHs = SpikeTrains_for_PSTHs;

    % store baseline/spontaneous ISIs (vector)
    rec_data.(cellID).ISI_baseline_vec = ISI_baseline_vec;

    % calculate and store baseline ISI probability density function (PDF)
    [ISI_pdf_x, ISI_pdf_y] = EmpiricalPDF(ISI_baseline_vec);
    rec_data.(cellID).ISI_pdf_y = ISI_pdf_y;
    rec_data.(cellID).ISI_pdf_x = ISI_pdf_x;
    % ... peak of the baseline ISI PDF ("ISI Peak" in Shin & Moore, 2019)
    [ISI_pdf_peak_y, ISI_pdf_peak_ind] = max(ISI_pdf_y);
    ISI_pdf_peak_x = ISI_pdf_x(ISI_pdf_peak_ind);
    rec_data.(cellID).ISI_pdf_peak_xy = [ISI_pdf_peak_x, ISI_pdf_peak_y];

    % ... coefficient of variation (CV) of baseline ISIs
    rec_data.(cellID).ISI_baseline_CV = std(ISI_baseline_vec) / mean(ISI_baseline_vec);

    % ... baseline Fano Factor
    % rec_data.(cellID).FanoFactor_baseline = var(FRs_baseline_vec) / mean(FRs_baseline_vec); % the wacky Shin way
    rec_data.(cellID).FanoFactor_baseline = var(double(spikeTrains_baseline_ms)) / mean(double(spikeTrains_baseline_ms),1);
    % note that the var(...) and mean(...) both operate down rows. What we
    % ultimately get is a mean of Fano Factors that were calculated per
    % millisecond.

    %% Create raw PSTHs as 2D matrices: trial_type x time (milliseconds)
    % ... and also Gaussian convolved PSTHs in the same format.
    % ... and also peak instantaneous evoked FRs in the user-defined time
    %     window (0 to user_time_cutoff_ms).
    % ... and also the average FR during baseline (1:onset-1) and
    %     stim (onset:onset+user_time_cutoff_ms), both from the raw PSTH.
    % ... and also all of the first spike latency stuff for each trial type
    %     (the probability density function, overall FSL, and FSL
    %     reliability).
    % ... and also average evoked FRs (binned) by averaging binned FRs
    %     across trials of each type.
    % ... and also Fano Factor based on spike counts per ms in the 
    %     user-defined analysis window, because an FR-based FF gave results 
    %     that were just bonkers.
    PSTHs_raw = [];
    PSTHs_conv = [];
    PeakEvokedFR = [];
    PeakEvokedFR_Latency = [];
    DeflectionResponseRatios = []; % [trial_type x deflection_num] matrix of individual deflection responses normalized to the first peak
    MeanFR_inst_baselines = []; % baselines of each trial type (vector), will be averaged to get MeanFR_inst_baseline (number)
    MeanFR_inst_stim = [];
    %FirstSpikeLatency_mean = [];
    MeanFR_stim = []; % separate means for each trial tag/type
    FanoFactor_stim = [];
    for trialTagInd = 1:length(uniqueTrialTags)
        SpikeTrains_trialType_ms = SpikeTrains_for_PSTHs{1,trialTagInd};
        SpikeTrains_trialType_stim_ms = SpikeTrains_trialType_ms(:, stim_duration_ms+1:stim_duration_ms+1+user_time_cutoff_ms);
        % ^ spike trains in milliseconds for the user-defined post-stim analysis window in this trial type

        PSTH_raw = mean(SpikeTrains_trialType_ms,1) * 1000; % *1000 to convert to Hz (spikes per second)
        PSTH_conv = conv(PSTH_raw, tempfilter, 'same');
        [PeakEvokedFR_ii, peak_ind] = max(PSTH_conv(stim_duration_ms+1:stim_duration_ms+1+user_time_cutoff_ms));
        % ^ note that stim_duration_ms+1 is the stim onset (inclusive), and
        % peak_ind is relative to that stim onset (= index 1).

        % Instantaneous FR means
        MeanFR_inst_baseline_ii = mean(PSTH_raw(1:stim_duration_ms));
        MeanFR_inst_stim_ii = mean(PSTH_raw(stim_duration_ms+1:stim_duration_ms+1+user_time_cutoff_ms));

        % Deflection response ratios
        if length(uniqueTrialTags) == 4 % if stimulation type is whisker
            DeflectionResponseRatios(trialTagInd,:) = GetDeflectionResponseRatios(PSTH_conv,relative_time_ms,MeanFR_inst_baseline_ii);
            %DeflectionResponseRatios(trialTagInd,:) = GetDeflectionResponseRatios(PSTH_raw,relative_time_ms,MeanFR_inst_baseline_ii);
        end

        % ... First Spike Latency PDF
        [FSL_pdf_x_ii, FSL_pdf_y_ii] = EmpiricalPDF(FirstSpikeLatency_perTrial{1,trialTagInd});
        FSL_pdf_x{1,trialTagInd} = FSL_pdf_x_ii;
        FSL_pdf_y{1,trialTagInd} = FSL_pdf_y_ii;
        % ... peak of the First Spike Latency PDF ("1st Spike Latency Reliability" in Shin
        %     & Moore, 2019) and its position in time ("1st Spike Latency"
        %     of the unit)
        [FSL_pdf_peak_y_ii, FSL_pdf_peak_ind] = max(FSL_pdf_y_ii);
        FSL_pdf_peak_x_ii = FSL_pdf_x_ii(FSL_pdf_peak_ind);
        FSL_pdf_peak_x(trialTagInd,1) = FSL_pdf_peak_x_ii; % unit's overall first spike latency for this trial type
        FSL_pdf_peak_y(trialTagInd,1) = FSL_pdf_peak_y_ii; % unit's first spike latency reliability for this trial type

        PSTHs_raw(trialTagInd,:) = PSTH_raw;
        PSTHs_conv(trialTagInd,:) = PSTH_conv;
        PeakEvokedFR(trialTagInd,1) = PeakEvokedFR_ii;
        PeakEvokedFR_Latency(trialTagInd,1) = peak_ind - 1; % remember, peak_ind is relative to the onset where onset_ind = 1
        MeanFR_inst_baselines(trialTagInd,1) = MeanFR_inst_baseline_ii;
        MeanFR_inst_stim(trialTagInd,1) = MeanFR_inst_stim_ii;
        %FirstSpikeLatency_mean(trialTagInd,1) = mean(FirstSpikeLatency_perTrial{1,trialTagInd}, 'omitnan');
        MeanFR_stim(trialTagInd,1) = mean(FRs_stim{1,trialTagInd});

        % Fano Factor (now using spike counts per millisecond instead of
        % firing rate)
        %FanoFactor_stim(trialTagInd,1) = var(FRs_stim{1,trialTagInd}) / mean(FRs_stim{1,trialTagInd}); % The wacky Shin method
        FanoFactor_stim(trialTagInd,1) = var(double(SpikeTrains_trialType_stim_ms)) / mean(double(SpikeTrains_trialType_stim_ms),1);
        % ^ this ends up being equivalent to the mean (ignoring NaN values)
        % of the elementwise division var(...) ./ mean(..., 1), and
        % accounts for variability not only in the *rate* of firing but
        % also in the precise timing of the spikes relative to stim onset.
    end
    rec_data.(cellID).PSTHs_raw = PSTHs_raw;
    rec_data.(cellID).PSTHs_conv = PSTHs_conv;
    rec_data.(cellID).PeakEvokedFR = PeakEvokedFR;
    rec_data.(cellID).PeakEvokedFR_Latency = PeakEvokedFR_Latency;
    rec_data.(cellID).MeanFR_inst_stim = MeanFR_inst_stim;
    rec_data.(cellID).MeanFR_inst_baseline = mean(MeanFR_inst_baselines);
    %rec_data.(cellID).FirstSpikeLatency_mean = FirstSpikeLatency_mean;
    rec_data.(cellID).MeanFR_stim = MeanFR_stim; % vector indexed by stim types/intensities, which are saved in trialTagsLabels.mat.
    rec_data.(cellID).FanoFactor_stim = FanoFactor_stim;

    rec_data.(cellID).FirstSpikeLatency_pdf_y = FSL_pdf_y;
    rec_data.(cellID).FirstSpikeLatency_pdf_x = FSL_pdf_x;
    rec_data.(cellID).FirstSpikeLatency = FSL_pdf_peak_x;
    rec_data.(cellID).FirstSpikeLatency_Reliability = FSL_pdf_peak_y;

    rec_data.(cellID).DeflectionResponseRatios = DeflectionResponseRatios;

    %% Calculate stimulus probability --> classify as responsive or non-responsive
    rec_data.(cellID) = ClassifyResponsivity_OneUnit(rec_data.(cellID));

end

end

