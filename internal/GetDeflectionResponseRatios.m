function DeflectionResponseRatios = GetDeflectionResponseRatios(PSTH_n,relative_time_ms,baselineFR)
%GETDEFLECTIONRESPONSERATIOS Summary of this function goes here
%   Detailed explanation goes here

%% Method 1: peaks with no baseline subtraction
% DeflectionResponseRatios = zeros(1,10);
% for deflectionNum = 1:10
%     start_ms = 50*(deflectionNum-1);
%     %end_ms = 50*deflectionNum;
%     end_ms = start_ms+30; % 30 ms window
% 
%     start_ind = find(relative_time_ms == start_ms);
%     end_ind = find(relative_time_ms == end_ms);
% 
%     DeflectionResponseRatios(deflectionNum) = max(PSTH_n(start_ind:end_ind));
% end
% 
% DeflectionResponseRatios = DeflectionResponseRatios / DeflectionResponseRatios(1); % normalization to first peak

%% Method 2: mean instantaneous FR per deflection bin with baseline subtraction
% DeflectionResponseRatios = zeros(1,10);
% for deflectionNum = 1:10
%     start_ms = 50*(deflectionNum-1);
%     %end_ms = 50*deflectionNum;
%     end_ms = start_ms+20; % doing 20 ms windows to be consistent with how the first peak was found
% 
%     start_ind = find(relative_time_ms == start_ms);
%     end_ind = find(relative_time_ms == end_ms);
% 
%     DeflectionResponseRatios(deflectionNum) = mean(PSTH_n(start_ind:end_ind));
% end
% 
% DeflectionResponseRatios = DeflectionResponseRatios - baselineFR; % baseline subtraction
% DeflectionResponseRatios = DeflectionResponseRatios / DeflectionResponseRatios(1); % normalization to first peak

%% Method 3: peaks, with subtraction of the minimum value in their deflection bin to (roughly) remove background tonic firing
% So the metric is the *magnitude* of the peak response relative to local
% background activity, which could include recurrence from previous
% deflections.
DeflectionResponses = zeros(1,10);
for deflectionNum = 1:10
    start_ms = 50*(deflectionNum-1);
    %end_ms = 50*deflectionNum;
    end_ms = start_ms+30; % 30 ms window

    start_ind = find(relative_time_ms == start_ms);
    end_ind = find(relative_time_ms == end_ms);

    DeflectionResponses(deflectionNum) = max(PSTH_n(start_ind:end_ind)) - min(PSTH_n(start_ind:end_ind));
end

DeflectionResponseRatios = DeflectionResponses / DeflectionResponses(1); % normalization to first peak magnitude

end

