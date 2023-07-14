function frequency = find_number_of_frequencies_fun(onset, offset, before_stimulation, after_stimulation, ...
            sampling_factor, frequency_count_identification_smoothening_method, analog_signal_optostim)
   
% Particular analog signal optostim
particular_analog_signal_optostim = ...
    analog_signal_optostim(onset - before_stimulation * sampling_factor:offset + after_stimulation * sampling_factor);
smoothened_analog_signal_for_frequency_identification = smoothdata(particular_analog_signal_optostim, frequency_count_identification_smoothening_method);

% Find the peaks in the data for a particular range which is from onset and offset from BlackRock
[peaks, peak_locations] = findpeaks(smoothened_analog_signal_for_frequency_identification);

threhold_peaks = 0.5e4; % uV
peak_locations_lower_than_threshold = find(peaks < threhold_peaks);
peak_locations(peak_locations_lower_than_threshold) = [];

% Count the number of peaks in the data
frequency = length(peak_locations)*2;

end
