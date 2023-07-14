function filtered_electrode_data = filter_data_fun(electrode_data_reordered, file_directory_to_save_data)


%% Butterworth bandpass filter design

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

filtered_electrode_data = filtfilt(B, A, electrode_data_reordered);

end

