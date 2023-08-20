function all_data = ReclassifyResponsivity(all_data, millisecondsAfterStim, alpha)
% For each unit, classify its sensory responsivity based on these
% parameters. Same method as in Shin & Moore, 2019.
%
% Assumes 30 kHz sampling rate.
%
% INPUTS
% - all_data
% - timeAfterStim: milliseconds after stimulation onset to use spike counts
%       for classification.
% - alpha: 0.05 or 0.01 for 95% or 99% CIs that the AUROC must exceed to be
%       counted as responsive.

groupNames = fieldnames(all_data);
for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};

    recNames = fieldnames(all_data.(groupName));
    for recNum = 1:length(recNames)
        recName = recNames{recNum};

        cellIDs = fieldnames(all_data.(groupName).(recName));
        for cellID_num = 1:length(cellIDs)
            cellID = cellIDs{cellID_num};
            
            stim_intensity = all_data.(groupName).(recName).(cellID).Stim_Intensity;

            % If there are multiple stimulation intensities, assume that the lowest is
            % 'catch' (zero) and the highest is 'max' & use 'catch' trials as baseline.
            % If there is only one stimulation intensity, use baseline FRs in one
            % column and stim FRs in the other.
            if length(unique(stim_intensity)) > 1
                inds_zero = find(stim_intensity == min(stim_intensity));
                inds_max = find(stim_intensity == max(stim_intensity));

                trialFRs_baseline = sum(all_data.(groupName).(recName).(cellID).SpikeTrains_stim(inds_zero,:), 2) / (size(all_data.(groupName).(recName).(cellID).SpikeTrains_stim(inds_zero,:),2) / 30000);
                trialFRs_stim = sum(all_data.(groupName).(recName).(cellID).SpikeTrains_stim(inds_max,1:millisecondsAfterStim*30), 2) / (size(all_data.(groupName).(recName).(cellID).SpikeTrains_stim(inds_max,1:millisecondsAfterStim*30),2) / 30000);
            else
                % trialFRs_baseline = unit_data.FRs_baseline';
                % trialFRs_stim = unit_data.FRs_stim';
                trialFRs_baseline = sum(all_data.(groupName).(recName).(cellID).SpikeTrains_baseline, 2)' / (size(all_data.(groupName).(recName).(cellID).SpikeTrains_baseline,2) / 30000);
                trialFRs_stim = sum(all_data.(groupName).(recName).(cellID).SpikeTrains_stim(:,1:millisecondsAfterStim*30), 2)' / (size(all_data.(groupName).(recName).(cellID).SpikeTrains_stim(:,1:millisecondsAfterStim*30),2) / 30000);
            end

            pred = [trialFRs_baseline; trialFRs_stim]; % predictor variable is FR

            baseline_labels = zeros(length(trialFRs_baseline),1);
            stim_labels = ones(length(trialFRs_stim),1);
            resp = [baseline_labels; stim_labels]; % response variable: 0 for baseline, 1 for stim

            %% fit logistic regression model and get the AUROC (w/ bootstrapped CI)
            warning('off','all');
            model = fitglm(pred,resp, 'Distribution','binomial', 'Link','logit');
            warning('on','all');
            
            scores = model.Fitted.Probability; % probability estimates
            
            % for 95% CI, pass 'Alpha',0.05 to perfcurve
            % for 99% CI, pass 'Alpha',0.01 to perfcurve
            [X,Y,T,AUC] = perfcurve(resp, scores, 1, 'NBoot',1000, 'BootType','bca', 'Alpha',alpha);
            
            
            % store Stimulus Probability (AUROC)
            all_data.(groupName).(recName).(cellID).StimProb = AUC; % [AUROC CI_lower CI_upper]

            %% calculate and store modulation index (MI): -1 <= MI <= 1
            % The MI tells us how strongly the unit is modulated by the stim, and in 
            % which direction. MI = (FR_stim - FR_baseline) / (FR_stim + FR_baseline)
            % Using mean FRs from the Gaussian-smoothed PSTH.
            FR_stim = mean(trialFRs_stim); % using evoked FR from the 'max' / last trial type
            FR_base = mean(trialFRs_baseline);
            ModulationIndex = (FR_stim - FR_base) / (FR_stim + FR_base);
            all_data.(groupName).(recName).(cellID).ModulationIndex = ModulationIndex;

            %% classify unit: -1 = negatively modulated, 0 = unresponsive, 1 = positively modulated
            % contra Shin & Moore, this has nothing to do with the sign of the AUROC!
            % AUROC only says that it's modulated in some direction; modulation index
            % determines which direction.
            if AUC(2) > 0.5 % unit is stimulus-responsive (i.e., evoked FR predicts stimulation condition)
                if ModulationIndex > 0 % unit is positively modulated
                    all_data.(groupName).(recName).(cellID).StimResponsivity = 1;
                else % MI is negative --> unit is negatively modulated
                    all_data.(groupName).(recName).(cellID).StimResponsivity = -1;
                end
            else % unit is not significantly stimulus-modulated/responsive
                all_data.(groupName).(recName).(cellID).StimResponsivity = 0;
            end
        end
    end
end

end

