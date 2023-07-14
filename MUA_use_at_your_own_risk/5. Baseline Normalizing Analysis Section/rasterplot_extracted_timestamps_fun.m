function rasterplot_extracted_timestamps_fun(extracted_time_stamps_NOT_normalized, electrodes_order, ...
    n_electrodes, n_trials, before_timestamp_range, after_timestamp_range, time_stamps, key_word, file_directory)

%{
    Plot the extracted timestamps figures separately to observe how the
    stimulation at different trials affect the neural activity of a mouse.
    These figures will be stored inside the sub-folder 
    "Extracted_Time_Stamps" within the folder "figures".
%}

new_electrode_number = n_electrodes;
original_electrode_number = electrodes_order(new_electrode_number);

% % Save all the figures to a separate folder called "figures/Extracted_Time_Stamps/".
% check_if_folder_exists(strcat(file_directory, 'figures/Extracted_Time_Stamps/', key_word, '/Electrode', string(original_electrode_number)));
    
% Plot the n_trials event
trial_number = n_trials;
stimulation_time = time_stamps(trial_number, 1);
stimulation_exposure_duration = time_stamps(1, 2)-time_stamps(1, 1);

left_time_bound = -before_timestamp_range+stimulation_time;
right_time_bound = stimulation_time+after_timestamp_range+stimulation_exposure_duration-1;
time_of_analysis = left_time_bound:1:right_time_bound;


extracted_time_stamps_NOT_normalized(trial_number, :, n_electrodes)

% auxilliary parameters 
Ntrials = size(row_stime_stamps);
T = 


plot(time_of_analysis, extracted_time_stamps_NOT_normalized(trial_number, :, n_electrodes), '-o');
xline(stimulation_time, '--r',{'Start of Stimulation'});
xline(stimulation_time+stimulation_exposure_duration, '--r',{'End of Stimulation'});
xlabel('Trial ' + string(n_trials) + '; Original Electrode ' + string(original_electrode_number))
ylabel('Number of spikes')

filename = strcat(file_directory, 'figures/Extracted_Time_Stamps/', key_word, '/Electrode', string(original_electrode_number), '/TrialNumber_', string(n_trials), '.png');
saveas(gcf, filename);

end

