% st_sim                        Spike timing analysis of time-domain WN analysis for simulated spike trains
%
% call                          [ stats ] = st_sim
% 15-Aug-22 AL
function [ stats ] = st_sim

% Load simulated spike trains
A = load('sim_spikes.mat');
sim_spikes = A.sim_spikes;
num_st = length(sim_spikes);
for i = 1:num_st
s = sim_spikes{i};
    Fs = 2500;
    % load the WN pattern
    WN                              = load( 'WN10k' );                          % frozen 1 s WN; mean 0, SD 1
    % resample and translate WN for stimulation:
    WN.template                     = WN.WN;
    wn0                             = resampleSig( WN.template, WN.Fs, Fs );

    % analyze
    [ statsi]              = stAnalysis( s, wn0 * ones( 1, size(s, 2 ) ) ...
        , 'ftype', 'wiener', 'graphics', [ 1 1  ], 'doRelLocal', 1, 'nreps1', 10, 'nreps2', 0, 'jmode', 1 );
    if i==1
        stats.HighP_rate = statsi;
    elseif i==2
        stats.HighP_temporal = statsi;
    elseif i==3
        stats.lowP_rate = statsi;
    elseif i==4
        stats.lowP_temporal = statsi;
    end

end
return