% st_neuronal                        Spike timing analysis of time-domain WN analysis for recorded spike trains
%
% call                          [ stats ] = st_neuronal
% 15-Aug-22 AL
function [ stats ] = st_neuronal

% Load recorded spike trains
A = load('neuronal_spikes.mat');
neuronal_spikes = A.neuronal_spikes;
num_st = length(neuronal_spikes);
for i = 1:num_st
s = neuronal_spikes{i};
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
        stats.DAPYR = statsi;
    elseif i==2
        stats.DAINT = statsi;
    elseif i==3
        stats.IDAINT = statsi;
    end

end
return