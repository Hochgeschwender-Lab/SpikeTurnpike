function [ccg_out,AB_sim_corrected,t] = JitterCorrectedCCG(spikeTimes1,spikeTimes2,duration_ms,halfbins,nreps,smoothingwin,plot_YorN)
% For two units' spike times, this function computes a jitter-corrected
% CCG using the interval jitter method (Platkiewicz, Stark, & Amarasingham,
% 2017, Neural Comput.) and detects/classifies whether there is a
% significant monosynaptic connection using simultaneous 95% acceptance
% bands, which correct for multiple comparisons. See Amarasingham et al., 
% 2011 (https://arxiv.org/pdf/1111.4296.pdf) for details about the latter.
%
% INPUTS
% - spikeTimes1, spikeTimes2: vectors of spike times in integer
%       milliseconds for two units (let's call them A and B)
% - duration_ms: duration of the time period or recording in milliseconds.
%       Used to determine the size of matrices for jittering. Defaults to
%       the latest spike time.
% - halfbins: milliseconds on either side of lag 0 to extend the CCG, so
%       length(ccg_out) == 2*halfbins + 1. Defaults to 100.
% - nreps: number of interval jitter CCGs to construct and average.
%       Defaults to 1000.
% - plot_YorN: 1 to plot the raw, mean jitter, and jitter-corrected CCGs.
%
% OUTPUTS
% - ccg_out: [milliseconds lag x 1] jitter-corrected CCG
% - t: [milliseconds lag x 1] time axis for ccg_out
% - sig: 1 if there is a significant excitatory A->B connection, -1 if
%       there is a significant inhibitory A->B connection, 0 otherwise.
% - AB_sim_corrected: simultaneous acceptance bands (col 1 is lower, col 2
%       is upper)
% - dcw: Directed Connection Weight as in Jia et al., 2021; the sum of the
%       jitter ccg in the interval (0,13) minus the sum of the jitter ccg
%       in the interval (-13,0).

if nargin < 7 || isempty(plot_YorN)
    plot_YorN = 0;
end
if nargin < 6
    smoothingwin = [];
end
if nargin < 5 || isempty(nreps)
    nreps = 1000;
end
if nargin < 4 || isempty(halfbins)
    halfbins = 100;
end
if nargin < 3 || isempty(duration_ms)
    duration_ms = max([max(spikeTimes1) max(spikeTimes2)]);
end

T = [spikeTimes1 spikeTimes2];
G = [repmat(1,[1 length(spikeTimes1)]) repmat(2,[1 length(spikeTimes2)])];

% normalization term: geometric mean of units' firing rates
% note: dividing by the millisecond duration to get Hz because the CCG bins
% are 1 ms
FRs = [length(spikeTimes1)/(duration_ms/1000) length(spikeTimes2)/(duration_ms/1000)];
norm = geomean(FRs);

% compute raw CCG
[ccg_raw, t, ~] = CCG(T, G, 1, halfbins, 1000, [1 2], 'count');
t = t';
% zero_lag_ind = find(t==0);
ccg_raw = ccg_raw(:,1,2) / norm; % only keep the CCG, not ACGs
if ~isempty(smoothingwin)
    ccg_raw = smoothdata(ccg_raw,1,"movmean",smoothingwin);
end

%% compute mean of nreps jittered CCGs
% first, jitter both spike trains as sparse matrices
x = zeros(duration_ms,2); % matrix of [time x units] for jittering
x(spikeTimes1,1) = 1;
x(spikeTimes2,2) = 1;
x = sparse(x);
x = stmix(x, 12, nreps, 2); % interval jitter with 25 ms width, nreps times

% generate nreps jitter CCGs and average them
jitterCCGs = zeros(length(ccg_raw),nreps);
for ii = 1:nreps
    jitteredST1 = find(x(:,ii));
    jitteredST2 = find(x(:,ii+nreps));
    
    T_jitter = [jitteredST1; jitteredST2]';
    G = [repmat(1,[1 length(jitteredST1)]) repmat(2,[1 length(jitteredST2)])];

    %T_jitter = [jitteredST1' spikeTimes2];
    %G = [repmat(1,[1 length(jitteredST1)]) repmat(2,[1 length(spikeTimes2)])];

    % Calculating a unique normalization term because the number of spikes
    % in jitteredST1 and jitteredST2 varies a bit, presumably because some
    % pairs of spikes are being jittered into the same millisecond
    % position. The effect of this is probably negligible, but normalizing
    % before averaging takes care of it anyway, just in case.
    FRs_ii = [length(jitteredST1)/(duration_ms/1000) length(jitteredST2)/(duration_ms/1000)];
    %FRs_ii = [length(jitteredST1)/(duration_ms/1000) FRs(2)];
    norm_ii = geomean(FRs_ii);

    jitterCCG_i = CCG(T_jitter, G, 1, halfbins, 1000, [1 2], 'count');
    jitterCCGs(:,ii) = jitterCCG_i(:,1,2) / norm_ii;
end
if ~isempty(smoothingwin)
    jitterCCGs = smoothdata(jitterCCGs,1,"movmean",smoothingwin);
end
meanJitterCCG = mean(jitterCCGs,2);

%% get jitter-derived pointwise 95% acceptance bands as described in https://arxiv.org/pdf/1111.4296.pdf
[AB_pt, AB_sim] = SurrogateAcceptanceBands(ccg_raw, jitterCCGs, 0.05);

% compute jitter-corrected CCG
ccg_out = ccg_raw - meanJitterCCG;
AB_pt_corrected = AB_pt - meanJitterCCG;
AB_sim_corrected = AB_sim - meanJitterCCG;

% %% define significant peak as in Siegle et al., 2021 (Nature), except
% % also capturing inhibitory connections.
% %noise_dist = ccg_out(t<=-50 | t>=50);
% %thres = [mean(noise_dist) + 7*std(noise_dist), mean(noise_dist) - 7*std(noise_dist)];
% [~,I_max] = max(abs(ccg_out(t>=-10 & t<=10))); % I = 11 is the same as t = 0
% t_max = I_max - 11;
% % [MinVal,I_min] = min(ccg_out(t>=-10 & t<=10));
% % t_min = I_min - 11;
% 
% %sig = (MaxVal > thres(1)) || (MinVal < thres(2)); % 1 if max exceeds upper threshold or min falls below lower threshold, 0 otherwise
% 
% % get the directed connection weight as in Jia et al., 2021
% dcw = sum(ccg_out(t>0 & t<13)) - sum(ccg_out(t>-13 & t<0)); % positive if A->B excitatory, negative if B->A excitatory, vice versa for inhibitory
% 
% % sig = 0;
% % if (MaxVal > thres(1)) && (I_max > zero_lag_ind) % significant excitatory A->B connection
% %     sig = 1;
% % elseif (MinVal < thres(2)) && (I_min > zero_lag_ind) % significant inhibitory A->B connection
% %     sig = -1;
% % end
% 
% % if (MaxVal > AB_sim_corrected(I_max,2)) && (t_max > 0) % significant excitatory A->B connection
% %     sig = 1;
% % elseif (MinVal < AB_sim_corrected(I_min,1)) && (t_min > 0) % significant inhibitory A->B connection
% %     sig = -1;
% % else
% %     sig = 0;
% % end
% 
% if (t_max > 1) && (ccg_out(t==t_max) > AB_sim_corrected(t==t_max,2))
%     sig = 1;
% elseif (t_max > 1) && (ccg_out(t==t_max) < AB_sim_corrected(t==t_max,1))
%     sig = -1;
% else
%     sig = 0;
% end


%% never mind, do it the cool way with simultaneous acceptance bands
% sig_peaks = find(ccg_out > AB_sim_corrected(:,2));
% sig_troughs = find(ccg_out < AB_sim_corrected(:,1));
% 
% if any(sig_peaks > zero_lag_ind & sig_peaks <= (zero_lag_ind+10)) % excitatory A->B
%     sig = 1;
% elseif any(sig_troughs > zero_lag_ind & sig_troughs <= (zero_lag_ind+10)) % inhibitory A->B
%     sig = -1;
% else
%     sig = 0;
% end

%% sanity check plots
if plot_YorN
    %noise_dist_raw = ccg_raw(t<=-50 | t>=50);
    %thres_raw = [mean(noise_dist_raw) + 7*std(noise_dist_raw), mean(noise_dist_raw) - 7*std(noise_dist_raw)];
    figure; T = tiledlayout(1,3);
    % nexttile(T); plot(t,ccg_raw); hold on; yline(thres_raw); xline([-10 10]); hold off; title(gca,'Raw CCG');
    % nexttile(T); plot(t,meanJitterCCG); hold on; xline([-10 10]); hold off; title(gca,'Mean Jitter CCG');
    % nexttile(T); plot(t,ccg_out); hold on; yline(thres); xline([-10 10]); hold off; title(gca,'Jitter-corrected CCG');
    nexttile(T); plot(t,ccg_raw,'k'); hold on; plot(t,AB_pt,'--y'); plot(t,AB_sim,'--r'); xline([-10 10]); hold off; title(gca,'Raw CCG');
    nexttile(T); plot(t,meanJitterCCG,'k'); hold on; plot(t,AB_pt,'--y'); plot(t,AB_sim,'--r'); xline([-10 10]); hold off; title(gca,'Mean Jitter CCG');
    nexttile(T); plot(t,ccg_out,'k'); hold on; plot(t,AB_pt_corrected,'--y'); plot(t,AB_sim_corrected,'--r'); xline([-10 10]); hold off; title(gca,'Jitter-corrected CCG');
    legend({'Data','Pointwise 95% Acceptance Band','','Simultaneous 95% Acceptance Band'}, 'Location','southoutside');
end

end