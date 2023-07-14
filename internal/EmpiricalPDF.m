function [pdf_x, pdf_y] = EmpiricalPDF(intervals_vec, bin_limits)
% Computes a sort of fake/numerical probability density function from a
% vector of either inter-spike intervals (one trial or concatenated 
% baselines) or first spike latencies (indexed by trial). The idea is to
% generate a "PDF" from undersampled data (e.g., units without enough 
% spikes to make a smooth histogram) and without specifying the 
% distribution a priori (e.g., MATLAB's pdf function).
% Uses histograms with 5 ms bins and 1 ms sliding step.
% Inferred protocol from Shin & Moore 2019, but they don't really describe
% it in detail.

%% Generate unnormalized histogram with 1 ms bins
if nargin < 2
    %bin_limits = [0 ceil(max(intervals_vec))];
    h = histogram(intervals_vec, 'BinWidth',1);
    %h = histogram(intervals_vec, 'BinMethod','integers', 'BinLimits',bin_limits);
else
    h = histogram(intervals_vec, 'BinWidth',1, 'BinLimits',bin_limits);
end

% Get x axis 'timeline' and y axis values
binEdges = h.BinEdges;
pdf_x = binEdges(1:end-1) + h.BinWidth/2;
y = h.Values;

%% Sum counts in 5 ms sliding window with 1 ms step
y_binned = movsum(y,5);

%% Normalize by bin size (5 ms) and n_observations to get a PDF
pdf_y = y_binned/(5*length(intervals_vec));

end

