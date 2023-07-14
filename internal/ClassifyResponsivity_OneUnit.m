function unit_data = ClassifyResponsivity_OneUnit(unit_data)
% Classify whether a unit is (1) positively responsive, (2) negatively
% responsive, or (3) not responsive to stimulation based on evoked firing
% rates.
% This is accomplished by ideal observer analysis (Shin & Moore, 2019) to
% determine whether a unit is sensory reponsive or non-responsive. Then,
% sensory responsive units are broken up by whether their mean change in
% firing rate is positive or negative.

stim_intensity = unit_data.Stim_Intensity;

% If there are multiple stimulation intensities, assume that the lowest is
% 'catch' (zero) and the highest is 'max' & use 'catch' trials as baseline.
% If there is only one stimulation intensity, use baseline FRs in one
% column and stim FRs in the other.
if length(unique(stim_intensity)) > 1
    % inds_catch = find(stim_intensity == min(stim_intensity));
    % inds_max = find(stim_intensity == max(stim_intensity));
    % trialFRs_baseline = unit_data.FRs_stim(inds_catch)';
    % trialFRs_stim = unit_data.FRs_stim(inds_max)';

    trialFRs_baseline = unit_data.FRs_stim{1,1};
    trialFRs_stim = unit_data.FRs_stim{1,end};
else
    % trialFRs_baseline = unit_data.FRs_baseline';
    % trialFRs_stim = unit_data.FRs_stim';
    trialFRs_baseline = unit_data.FRs_baseline_vec;
    trialFRs_stim = unit_data.FRs_stim{1,1};
end

pred = [trialFRs_baseline; trialFRs_stim]; % predictor variable is FR

baseline_labels = zeros(length(trialFRs_baseline),1);
stim_labels = ones(length(trialFRs_stim),1);
resp = [baseline_labels; stim_labels]; % response variable: 0 for baseline, 1 for stim

%% fit logistic regression model and get the AUROC (w/ bootstrapped CI)
model = fitglm(pred,resp, 'Distribution','binomial', 'Link','logit');

scores = model.Fitted.Probability; % probability estimates

% for 95% CI, pass 'Alpha',0.05 to perfcurve
% for 99% CI, pass 'Alpha',0.01 to perfcurve
[X,Y,T,AUC] = perfcurve(resp, scores, 1, 'NBoot',1000, 'BootType','bca', 'Alpha',0.01);


% store Stimulus Probability (AUROC)
unit_data.StimProb = AUC; % [AUROC CI_lower CI_upper]

%% calculate and store modulation index (MI): -1 <= MI <= 1
% The MI tells us how strongly the unit is modulated by the stim, and in 
% which direction. MI = (FR_stim - FR_baseline) / (FR_stim + FR_baseline)
% Using mean FRs from the Gaussian-smoothed PSTH.
FR_stim = unit_data.MeanFR_inst_stim(end); % using evoked FR from the 'max' / last trial type
FR_base = unit_data.MeanFR_inst_baseline;
ModulationIndex = (FR_stim - FR_base) / (FR_stim + FR_base);
unit_data.ModulationIndex = ModulationIndex;

%% classify unit: -1 = negatively modulated, 0 = unresponsive, 1 = positively modulated
% contra Shin & Moore, this has nothing to do with the sign of the AUROC!
% AUROC only says that it's modulated in some direction; modulation index
% determines which direction.
if AUC(2) > 0.5 % unit is stimulus-responsive (i.e., evoked FR predicts stimulation condition)
    if ModulationIndex > 0 % unit is positively modulated
        unit_data.StimResponsivity = 1;
    else % MI is negative --> unit is negatively modulated
        unit_data.StimResponsivity = -1;
    end
else % unit is not significantly stimulus-modulated/responsive
    unit_data.StimResponsivity = 0;
end

%% sanity check plot of ROC
% figure
% plot(X,Y)
% xlabel('False positive rate') 
% ylabel('True positive rate')
% title('ROC for Classification by Logistic Regression')
% 
% AUC

end

