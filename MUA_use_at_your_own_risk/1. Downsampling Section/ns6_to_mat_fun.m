function electrode_data_original = ns6_to_mat_fun(ns6_file_directory)
%{
    This file is for getting a blackrock ns6 to an int16 matrix of 
    [nChannels x nEvents] and saving it as a .mat for further analysis.
%}

%% Read and save data
%{
 This function reads the data from the directory as a struct data type and
 returns the actual data back.
%}

% Read the .ns6 data
ns6_struct = openNSx('read','c:1:32', ns6_file_directory); %only will read channels 1-32, 33,34,35 are the analog signals

% Return the data from within the struct data
electrode_data_original = ns6_struct.Data;

end

