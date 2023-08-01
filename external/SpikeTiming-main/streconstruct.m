% streconstruct         reconstruct analog waveform from spike trains
% 
% call                  [ yhat, stats, xc, ax ] = streconstruct( x, y )
%
% gets                  x               spike train/s. sparse matrix (m x n) or a vector of times
%                       y               analog signal. m x 1 (or m x n, or p x 1)
%
% returns               yhat            reconstructed analog signal, same dimensions as y
%                       f               reconstruction filters (mf x nf matrix)
%                       xc              classified spike train/s (same dimension as x)
%                       ax              axes handles (on a new figure)
%
% optional arguments (given as name/value pairs)
%   
%                       minSpikes             {5}; the minimal total number of spikes in x
%                       minTrials             {4}; the minimal number of trials
%                       minSpikesPerTrial     {1}; the minimal mean number of spikes/trial
%                       prepy                 {'zs'} (standradize) or 'mean' (remove mean)
%                       prepMode              {'flat'} (ignore spikes at edges) or 'pad' (zero-pad; triangular-bias in STA/auto-corr)
%                       ftype                 first order:      {'sta'}, 'wiener', 'wsta', 'corr'
%                                             second order:     'csta' 'cwiener', 'cwsta'
%                                                               'xsta' 'xwiener', 'xwsta' (requires x2)
%                       x2                    {[]}; a reference set of spike trains (same arrangement as x). relevant for across-train, second order models only
%                       nfilt                 {50}; [samples]; filter range. a scalar indicates [ -1 1 ] * nfilt
%                       maxISIbins            {10}; maximal number of ISI bins (relevant for 'csta' and 'xsta' only)
%                       minCountPerBin        {10}; minimal number of spikes/ISI bin (relevant for 'csta' and 'xsta' only)
%                       isiedges              {[]}; [samples]; externaly-specified set of ISI bin edges (overrides maxISIbins/minCountPerBin)
%                       lagsr                 {0}; [samples]; range of lags to optimize over. a scalar
%                                                 indicates a single lag, and positive lags indicate causal filter
%                       cctype                {'rank'} or 'pmoment'; CC type used for post-hoc quantification 
%                       nreps                 {10} number of randomizations per trial
%                       rhalfwin              {[]} (shuffle) or a number (jitter half win)
%                       graphics              {1} or 0; to plot or not (one summary figure)
%                       Fs                    {1}; [samples/s]; sampling rate - used ONLY for ftim/graphics (i.e. if is kHz, abscissa in ms etc)
% 
% Notes:
% 1. This is a low level routine, so everything is in samples; x and y must be of the same sampling rate. 
%       If that is not the case, the calling routine should downsample/upsample appropriately
% 2. Only relevant for zero-mean signals (waveform reconstruction)
% 
% calls                 ParseArgPairs
%                       calc_fwhh, calc_pearson, calc_spearman
%                       isisort, my_xcorr, myjet, plot_raster
%                       sta, stacalc, stahat, stmix, stwiener
%                       zsc (mex)
%
% references            Bialek et al, 1991, Science
%                       Rieke et al., 1997, Spikes
%                       van Rossum 2001, Neural Computation
%                       Schreiber et al., 2003, Neurocomputing
%
% see also              streconstructXval

% 09-jun-14 ES

% last update
% August 2022

function [ yhat, stats, xc, ax ] = streconstruct( x, y, varargin )

yhat                            = [];
stats                           = [];
xc                              = [];
ax                              = [];

%----------------------------------------------------------------------
% constants
%----------------------------------------------------------------------
% reliability
BLOCKSIZE                       = 100;

% coherence
fROI                            = [ 0 200 ];
M                               = 0.25;
mtNW                            = 3;
dflag                           = '';

% graphics
colors                          = [ 0 0 0.7; 0 0 0; 1 0 0 ];                % x, y, yhat
blackColor                      = [ 0 0 0 ];
whiteColor                      = [ 1 1 1 ];
grayColor                       = [ 1 1 1 ] * 0.6;
purpleColor                     = [ 1 0.5 1 ]; 
redColor                        = [ 1 0 0 ];
tstr                            = '';

%----------------------------------------------------------------------
% arguments
%----------------------------------------------------------------------
nargs                           = nargin;
if nargs < 2 || isempty( x ) || isempty( y )
    return
end
[ minSpikes, minTrials, minSpikesPerTrial, prepy ...
    , prepMode ...
    , ftype, x2, nfilt ...
    , maxISIbins, minCountPerBin, isiedges ...
    , lagsr ...
    , cctype, nreps, rhalfwin, doRel ...
    , doScatter, doCoherence, Fs, graphics, verbose ] = ParseArgPairs(...
    { 'minSpikes', 'minTrials', 'minSpikesPerTrial', 'prepy' ...
    , 'prepMode' ...
    , 'ftype', 'x2', 'nfilt' ...
    , 'maxISIbins', 'minCountPerBin', 'isiedges' ...
    , 'lagsr' ...
    , 'cctype', 'nreps', 'rhalfwin', 'doRel' ...
    , 'doScatter', 'doCoherence', 'Fs', 'graphics', 'verbose' }...
    , { 5, 4, 1, 'zs' ...
    , 'flat' ...
    , 'sta', [], 50 ...
    , 10, 5, [] ...
    , 0 ...
    , 'rank', 10, [], 0 ...
    , 0, 0, 1, 0, 1 }...
    , varargin{ : } );

% scale analog
if ~ismatrix( y )
    y                           = y( : );
end
y0                              = y;
prepy                           = lower( prepy );
switch prepy
    case 'zs'
        y                       = zsc( y0, 1 );
        zstr                    = 'Z';
    case 'mean'
        y                       = bsxfun( @minus, y, nanmean( y, 1 ) );
        zstr                    = 'zraw';
    otherwise
        error( 'written ONLY for zero-mean signals' )
end

% input size
y1                              = y;
if issparse( x )
    [ m, nt ]                   = size( x );
    [ my, ny ]                  = size( y );
    if nt > 1
        if ny == 1
            if m == my
                y               = repmat( y, [ 1 nt ] );
            else
                error( 'input size mismatch: X sparse, Y mismatching vector!' )
            end
        elseif nt ~= ny
            error( 'input size mismatch: X sparse, Y matrix with mismatching columns!' )
        elseif m ~= my
            error( 'input size mismatch: X sparse, Y matrix with mismatching rows!' )
        end
    elseif m ~= my
        error( 'input size mismatch: X sparse, Y matrix with mismatching rows!' )
    end
else
    error( 'not supported yet: X must be sparse' )
end
if nt < minTrials
    if verbose
        fprintf( '%s: Too few trials!! (%d)\n', upper( mfilename ), nt )
    end
    return
end

% filter type
ftype                           = lower( ftype );
if ~ismember( ftype, { 'wiener', 'cwiener', 'xwiener' ...
        , 'sta', 'csta', 'xsta' ...
        , 'wsta', 'cwsta', 'xwsta' ...
        , 'corr' ...
        , 'sta1', 'sta2', 'sta3' } )
    error( 'mismatching filter type' )
end
switch ftype
    case { 'sta1' }
        ftype                   = 'sta';
    case { 'sta2' }
        ftype                   = 'csta';
    case { 'sta3' }
        ftype                   = 'xsta';
end
if ismember( ftype, { 'xsta', 'xwiener', 'xwsta' } ) && ~isequal( [ m nt ], size( x2 ) )
    error( 'input size mismatch: X and X2 must have identical dimensions' )
end

% filter length
if length( nfilt ) == 1
    nfilt                       = nfilt * [ -1 1 ];
end
if length( nfilt ) ~= 2
    error( 'NFILT must be a range' )
end
nfilt                           = round( sort( nfilt ) );               	% to be computed, may be a-causal
maxlag                          = max( abs( nfilt ) );
nfiltf                          = max( abs( nfilt ) ) * [ -1 1 ];       	% full, a-causal

% isis bins
if ismember( ftype, { 'cwiener', 'csta', 'cwsta' ...
        , 'xwiener', 'xsta', 'xwsta' } ) 
    if isempty( isiedges )
        if ( isempty( maxISIbins ) || maxISIbins < 1 )
            error( 'MAXISIBINS must be non-negative' )
        end
        maxISIbins              = round( maxISIbins );
    else
        [ ~, tmp ]              = size( isiedges );
        if tmp ~= 1
            error( 'ISIBINS must be a vector of edges' )
        end
        isiedges                = sort( isiedges );
    end
end

% randomization scale
rhalfwin                        = abs( round( rhalfwin ) );

%----------------------------------------------------------------------
% prepare for filter computation
%----------------------------------------------------------------------
% pad to prevent edge effects
switch prepMode
    case 'pad'
        padbuffer               = zeros( maxlag, nt );
        ypad                    = [ y; padbuffer ];                        	% does not remove any spikes
        xpad                    = [ x; padbuffer ];                      	% does not modify the analog
    case { 'flat', 'unbiased' }
        if nfilt( 1 ) < 0
            x( 1 : abs( nfilt( 1 ) ), : ) = 0;
        end
        if nfilt( 2 ) > 0
            mx                  = size( x, 1 );
            zi                  = ( mx - nfilt( 2 ) + 1 ) : mx;
            x( zi, : )          = 0;
        end
        xpad                    = x;                                     	% removed data from the spikes
        ypad                    = y;                                     	% no padding of the analog
        padbuffer               = [];                                    	% to support c-sta/c-wiener etc
end

% initialize
nspikes                         = full( sum( x( : ) ) );                    % number of spikes
nf                              = diff( nfilt ) + 1;
f                               = NaN * ones( nf, 1 );
fsem                            = f;
ftim                            = ( nfiltf( 1 ) : nfiltf( 2 ) )' / Fs;   	% [s]

bedges                          = [];
cpb                             = [];
sidx                            = sparse( nspikes, 1 );

% make sure enough spikes:
if nspikes < minSpikes || nspikes / nt < minSpikesPerTrial
    if verbose
        fprintf( '%s: Too few spikes!! (%d trials, %d spikes)\n', upper( mfilename ), nt, nspikes )
    end
    return
end

%----------------------------------------------------------------------
% compute filters
%----------------------------------------------------------------------

xc                              = x;
switch ftype
    
    case 'corr'                                                             % reverse correlation filter
        f                       = my_xcorr( ypad( : ), full( xpad( : ) ), maxlag, 0 );
        f                       = f / nspikes;
        fsem                    = ones( size( fsem ) ) / sqrt( nspikes );
        
    case 'sta'                                                              % STA filter (mathematically identical, just faster for sparse trains)
        xt                      = find( xpad( : ) ); 
        [ f, ss, nn ]           = sta( xt( : ), ypad( : ), nfilt );
        fsem                    = ss / sqrt( nn );

    case 'wsta'                                                             % input-whitened STA filter (first-order Wiener filter)
        xt                      = find( xpad( : ) );
        [ f, fsem ]             = wsta( xt, ypad( : ), nfilt );

    case 'wiener'                                                           % output-whitened 'Wiener' filter (point process implementation)
        xt                      = find( xpad( : ) );
        [ f, fsem ]             = stwiener( xt, ypad( : ), nfilt );
        
    case { 'cwiener', 'csta', 'cwsta' ...
            'xwiener', 'xsta', 'xwsta' }                                    % conditional w/STA/Wiener filter (depending on the preceding ISI)
        
        % sort the spikes by the ISIs:
        nspikesEff              = sum( max( full( sum( x, 1 ) ), 1 ) - 1 );
        if isempty( isiedges )
            nibins      = min( maxISIbins, floor( nspikesEff / minCountPerBin ) );
            if nibins < 2 
                if verbose
                    fprintf( '%s: Too few ISIs!! (%d trials, %d valid ISIs)\n' ...
                        , upper( mfilename ), nt, nspikesEff )
                end
                return
            end
            isiedges            = nibins;
        end
        switch ftype
            case { 'csta', 'cwiener', 'cwsta' }
                [ xc, bedges, cpb ] = isisort( x, isiedges );
            case { 'xsta', 'xwiener', 'xwsta' }
                [ xc, bedges, cpb ] = isisort( x, isiedges, x2 );
        end
        nibins                  = length( bedges ) - 1;
        sidx                    = full( xc( xc > 0 ) );
        
        % compute filter bank:
        xt                      = find( [ xc; padbuffer ] );
        switch ftype
            case { 'csta', 'xsta' }
                cmethod         = 'sta';
            case { 'cwsta', 'xwsta' }
                cmethod         = 'wsta';
            case { 'cwiener', 'xwiener' }
                cmethod         = 'wiener';
        end
        [ f, fsem ]             = stacalc( xt, ypad( : ), nfilt, sidx, 1 : nibins, cmethod );
        
end

% compute the 1st order STA filter anyhow
if ismember( ftype, { 'corr', 'sta', 'wiener', 'wsta' } )
    f1                          = f;
    f1sem                       = fsem;
else
    [ f1, ss, nn ]              = sta( xt( : ), ypad( : ), nfilt );
    f1sem                       = ss / sqrt( nn );
end

%----------------------------------------------------------------------
% reconstruction
%----------------------------------------------------------------------

% pad the filter to allow simple firfilt-based reconstruction
f0                              = f;
f0s                             = fsem;
f01                             = f1;
f01s                            = f1sem;
[ ~, fi ]                       = intersect( nfilt( 1 ) : nfilt( 2 ), nfiltf( 1 ) : nfiltf( 2 ) );
f                               = zeros( 2 * maxlag + 1, size( f, 2 ) );
f( fi, : )                      = f0;
fsem                            = zeros( 2 * maxlag + 1, size( f, 2 ) );
fsem( fi, : )                   = f0s;
f1                              = zeros( 2 * maxlag + 1, 1 );
f1( fi, : )                     = f01;
f1sem                           = zeros( 2 * maxlag + 1, 1 );
f1sem( fi, : )                  = f01s;

% determine the lag
if length( lagsr ) == 1                                                     % fixed lag
    optlag                      = lagsr;
    f                           = local_shift( f, -optlag, 0 );
    if optlag > 0
        f( ftim > 0, : )        = 0;
    end        
else
    % optimize over all possible lags:
    fprintf( '%s: Optimizing STA lag: ', upper( mfilename ) )
    lags                        = lagsr( 1 ) : lagsr( 2 );
    nlags                       = length( lags );
    fx                          = ftim .* Fs;
    fidx                        = fx <= 0;
    cci                         = zeros( 2 * maxlag + 1, nlags );
    for i                       = 1 : nlags
        fprintf( '.' )
        lag                     = lags( i );
        fhat                    = local_shift( f, -lag, 0 );                % shift
        fhat( ~fidx, : )        = 0;                                      	% make causal
        yhat                    = stahat( xc, fhat, size( y, 1 ), 0 );
        cci( :, i )             = mean( my_xcorr( y, yhat, maxlag, -1 ), 2 );
    end
    % keep the causal filter with the best performance:
    [ ~, lagidx ]               = find( cci == max( cci( : ) ) );
    optlag                      = lags( lagidx );
    f                           = local_shift( f, -optlag, 0 );
    f( ~fidx, : )               = 0;
    fprintf( 'Optimal lag: %d\n', optlag )
end

% reconstruct + shift back
yhat                            = stahat( xc, f, size( y, 1 ), optlag );

%----------------------------------------------------------------------
% quantification + statistics
%----------------------------------------------------------------------

% initialize
pvals                           = NaN * ones( 1, nt );
pvalAll                         = [];
rAll                            = [];

% reconstruction quality (Q := ccAll)
switch cctype
    case 'rank'
        yr                      = local_rankcols( y );
        [ cct, ftc ]            = my_xcorr( yr, local_rankcols( yhat ), maxlag, -1 );
        act                     = my_xcorr( yr( :, 1 ), yr( :, 1 ), maxlag, -1 );
        cctAll                  = my_xcorr( yr( : ), local_rankcols( yhat( : ) ), maxlag, -1 );
        ccstr                   = 'Rank CC';
    case 'pmoment'
        [ cct, ftc ]            = my_xcorr( y, yhat, maxlag, -1 );
        act                     = my_xcorr( y( :, 1 ), y( :, 1 ), maxlag, -1 );
        cctAll                  = my_xcorr( y( : ), yhat( : ), maxlag, -1 );
        ccstr                   = 'CC';
end
cc                              = cct( maxlag + 1, : );                     % zero lag cc
ccAll                           = nanmean( cc );                         	% mean over trials
[ ccmax, ccmaxlag ]             = max( cct, [], 1 );                     	% max cc - may be noisy
[ ccmaxAll, ccmaxAllLag ]       = max( cctAll, [], 1 );

% randomize the zero-lag cc (Q) for p-value
if nreps > 0
    fprintf( '%s: Randomizing (%d trials, %d times) ', upper( mfilename ), nt, nreps )
end
yhatRand                        = zeros( size( yhat, 1 ), nt );
for i                           = 1 : nt
    if nreps > 0
        fprintf( '.' )
    else
        continue
    end
    xm                          = stmix( x( :, i ), rhalfwin, nreps );
    switch ftype
        case { 'corr', 'sta', 'wiener', 'wsta' }                            % 1st order
            yhatR               = local_firfilt( full( xm ), f );
        case { 'cwiener', 'xwiener', 'csta', 'xsta', 'cwsta', 'xwsta' }     % 2nd order
            xm                  = isisort( xm, bedges, x2 );
            yhatR               = stahat( xm, f, size( y, 1 ), optlag );
    end
    yhatRand( :, i )            = yhatR( :, 1 );
    switch cctype
        case 'rank'
            r                   = calc_pearson( yr( :, i ) * ones( 1, nreps ), local_rankcols( yhatR ) );
        case 'pmoment'
            r                   = calc_pearson( y( :, i ) * ones( 1, nreps ), yhatR );
    end
    rAll                        = [ rAll r ];
    if nreps < 100                                                          % assume Gaussian distribution of the shuffled data
        mm                      = nanmean( r );
        ss                      = nanstd( r );
        pvals( i )              = 1 - normcdf( cc( i ), mm, ss );
    else                                                                    % compute empirical p-value
        pvals( i )              = ( sum( r >= cc( i ) ) + 1 ) / ( nreps + 1 );
    end
end
if nreps > 0 && nreps < 100
    mm                          = nanmean( rAll );
    ss                          = nanstd( rAll );
    pvalAll                     = 1 - normcdf( ccAll, mm, ss );
elseif nreps > 100
    pvalAll                     = ( sum( rAll >= ccAll ) + 1 ) / ( nreps + 1 );
end
if nreps > 0
    fprintf( '\n' )
end
if isempty( rAll )
    rAll                        = NaN; 
    pvalAll                     = NaN;
end 

% temporal statistics (T := lagSta)
[ widSta, lagSta ]              = calc_fwhh( f1 );                       	% from the 1st order filter
lagSta                          = lagSta - maxlag - 1;
wid                             = calc_fwhh( cct );                        	% from the reconstruction (one per trial)
widAll                          = calc_fwhh( cctAll );                    	% from the reconstruction (concatenated data)
widAct                          = calc_fwhh( act );                        	% upper bound: the FWHH of the input auto-correlation

% reliability + width by pair-wise CC between reconstructions (R := ccRel); Schreiber et al., 2003
widCC2                          = NaN;
ccRel                           = [ NaN NaN NaN ];
ccRelRand                       = [ NaN NaN NaN ];
if doRel
    
    aidx                        = local_make_ai( nt );
    blocks                      = local_makeblocks( size( aidx, 1 ), BLOCKSIZE );
    nblocks                     = size( blocks, 1 );
    cc2                         = zeros( 2 * maxlag + 1, 1 );
    for i                       = 1 : nblocks
        bidx                    = blocks( i, 1 ) : blocks( i, 2 );
        cc2i                    = my_xcorr( yhat( :, aidx( bidx, 1 ) ), yhat( :, aidx( bidx, 2  ) ), maxlag, -1 );
        cc2                     = cc2 + nanmean( cc2i, 2 );
    end
    cc2                         = cc2 / nblocks;
    widCC2                      = calc_fwhh( mean( cc2, 2 ) );              % width of the reliability cross-correlation
    rel                         = cc2( maxlag + 1, : );
    relPval                     = signrank( rel, 0 );
    if nanmedian( rel ) > 0
        relPval                 = relPval / 2;
    else
        relPval                 = 1 - relPval;
    end
    ccRel                       = [ nanmean( rel ) std( rel ) / sqrt( size( rel, 1 ) ) relPval ];
    
    if nreps > 0
        cc2rand                 = zeros( 1 );
        for i                   = 1 : nblocks
            bidx                = blocks( i, 1 ) : blocks( i, 2 );
            switch cctype
                case 'rank'
                    cc2i        = calc_spearman( yhat( :, aidx( bidx, 1 ) ), yhat( :, aidx( bidx, 2  ) ) );
                case 'pmoment'
                    cc2i        = calc_pearson( yhat( :, aidx( bidx, 1 ) ), yhat( :, aidx( bidx, 2  ) ) );
            end
            cc2rand             = cc2rand + nansum( cc2i, 2 );
        end
        cc2rand                 = cc2rand / nblocks;
        relRand                 = cc2rand;
        relRandPval             = signrank( relRand, 0 );
        if nanmedian( relRand ) > 0
            relRandPval         = relRandPval / 2;
        else
            relRandPval         = 1 - relRandPval;
        end
        ccRelRand               = [ nanmean( relRand ) ...
            std( relRand ) / sqrt( size( relRand, 1 ) ) relRandPval ];
    end
    
end

% rescale (for scatter, for plotting)
if isequal( prepy, 'zs' )
    yhat                        = zsc( yhat, 1 );
end

% scatter - plot <y>|yhat in [a,b) vs. yhat in [a,b); Rieke et al., 1997
if doScatter || graphics
    nybins                      = min( floor( m / 20 ), 50 );
    [ ~, ~, yidx ]              = local_make_equal_bins( yhat( : ), nybins );
    ysct                        = zeros( nybins, 2 );
    for i                       = 1 : nybins
        idx                     = yidx == i;
        ysct( i, : )            = [ mean( yhat( idx ) ) mean( y( idx ) ) ];
    end
    switch cctype
        case 'rank'
            ccLinear            = calc_spearman( ysct );
        case 'pmoment'
            ccLinear            = calc_pearson( ysct );
    end
else
    ccLinear                    = NaN;
end

% coherhence and phase between y and yhat
if doCoherence || graphics
    nFFT                        = 2^floor( log2( M * Fs * 1000 ) ); 
    nWindow                     = nFFT;
    [ yo, frq ]                 = local_mtcsd1( [ y( :, 1 ) yhat ], nFFT ...
        , Fs * 1000, nWindow, nWindow/2, mtNW, dflag );
    frqidx                      = frq >= min( fROI ) & frq <= max( fROI );
    yo                          = yo( frqidx, :, : );
    frq                         = frq( frqidx );
    p1                          = repmat( yo( :, 1, 1 ), [ 1 nt ] );
    p2                          = squeeze( yo( :, 1, 2 : ( nt + 1 ) ) );
    c12                         = squeeze( yo( :, 2, 2 : ( nt + 1 ) ) );
    cohs                        = single( abs( c12 .^ 2 ) ./ ( p1 .* p2 ) );
    phs                         = single( atan2( imag( c12 ), real( c12 ) ) );
else
    frq                         = [];
    cohs                        = [];
    phs                         = [];
end

%----------------------------------------------------------------------
% summarize results
%----------------------------------------------------------------------
% input:
stats.nspikes                   = nspikes;                                  % [count] number of spikes
stats.nt                        = nt;                                       % [count] number of trials
stats.Fs                        = Fs;                                       % [samples/s]
stats.rhalfwin                  = rhalfwin;                                 % [samples]

% filters:
stats.f                         = f;                                        % [Z] actual (STA/cSTA/xSTA) filter
stats.fsem                      = fsem;                                     % [Z] SEM of actual (1st or 2nd order) filter
stats.ftim                      = ftim;                                     % [s] for filters
stats.f1                        = f1;                                       % [Z] 1st order (STA) filter
stats.f1sem                     = f1sem;                                    % [Z] SEM of 1st order filter
stats.isiedges                  = bedges;                                   % [samples] ISI bin edges
stats.cpb                       = cpb;                                      % [count] per ISI bin
stats.sidx                      = sidx;                                     % [class] ISI class of each spike
stats.optlag                    = optlag;                                   % [samples] optimized (or fixed) lag (how much the filter+reconstruction was shifted) 

% reconstruction quality:
stats.ccLinear                  = ccLinear;                                 % quantifies how bad are saturation effects (linear regression of yhat vs. y)
stats.frq                       = frq;                                      % [Hz]; frequency bins for the coherence/phase estimates
stats.cohs                      = cohs;                                     % nfbins x nt; coherence between analog and its reconstruction
stats.phs                       = phs;                                      % [rad]

% summary stats
stats.cc                        = cc;                                       % [cc] @ zero-lag, 1 x nt
stats.pvals                     = pvals;                                    % per trial (by randominzation - {shuffle} or jitter)
stats.ccAll                     = [ ccAll std( cc ) / sqrt( size( cc, 1 ) ) ];                % [cc] @ zero-lag, mean and SEM over trials
stats.pval                      = pvalAll;                                  % same, for entire dataset
stats.ccRand                    = [ nanmean( rAll ) std( rAll ) / sqrt( size( rAll, 1 ) ) ]; % randomized

stats.ccRel                     = ccRel;                                    % zero-lag CC between reconstructions, mean and SEM
stats.ccRelRand                 = ccRelRand;                                % randomized

stats.lagsta                    = lagSta;                                   % [samples], here negative are causal
stats.widsta                    = widSta;                                   % [samples]; from the filter
stats.wid                       = wid;                                      % [samples] FWHH: from the reconstruction (for each trial)
stats.widAll                    = widAll;                                   % [samples] FWHH: from the reconstruction (concatenated data)
stats.widCC2                    = widCC2;                                   % [samples] FWHH: from the CC between reconstructions
stats.widAct                    = widAct;                                   % [samples] from the input

stats.ccmax                     = ccmax;                                    % [cc] either rank or pearson (for each trial)
stats.ccmaxlag                  = ccmaxlag - maxlag;                        % [samples]; lag of max cc (for each trial)
stats.ccmaxAll                  = ccmaxAll;                
stats.ccmaxAllLag               = ccmaxAllLag - maxlag;

% the "most important" results are:
% 1.1. ccAll        Q:  CC (trial-averaged mean and SEM) at zero-lag, for the selected filter (STA/cSTA etc)
% 1.2. pval             pval of the ccAll
% 2. lagsta         T:  [samples] time lag, based on the STA filter
% 3. ccRel          R:  CC (trial-averaged mean, SEM, and p-value) for the selected filter
% 4. halfwid        [samples] half width, based reconstructed concatenated data + input auto-correlation

%----------------------------------------------------------------------
% summarize graphically
%----------------------------------------------------------------------
if graphics 
    ytim                        = ( 1 : m )' / Fs;
    
    % layout
    figure
    ax( 1 )                     = axes( 'position', [ 0.05  0.75 0.9 0.2 ] );       % spike trains
    ax( 2 )                     = axes( 'position', [ 0.05  0.55 0.9 0.2 ] );       % analog + yhat
    ax( 3 )                     = axes( 'position', [ 0.05  0.3  0.2 0.2 ] );       % STA
    ax( 4 )                     = axes( 'position', [ 0.3   0.3  0.2 0.2 ] );       % (cSTA)
    ax( 5 )                     = axes( 'position', [ 0.525 0.3  0.2 0.2 ] );       % y-yhat xcc
    ax( 6 )                     = axes( 'position', [ 0.75  0.3  0.2 0.2 ] );       % y-yhat scatter
    ax( 7 )                     = axes( 'position', [ 0.05  0.05 0.2 0.2 ] );       % 2D optimization
    ax( 8 )                     = axes( 'position', [ 0.3   0.05 0.2 0.2 ] );       % 1D optimization
    ax( 9 )                     = axes( 'position', [ 0.525 0.05 0.2 0.2 ] );       % coherence
    ax( 10 )                    = axes( 'position', [ 0.75 0.05 0.2 0.2 ] );        % pvals
    
    % output - spikes
    subplot( ax( 1 ) )
    plot_raster( x, ytim, [], [], [], colors( 1, : ) );

    % reconstruction example
    subplot( ax( 2 ) )
    [ maxval, maxidx ]          = max( cc );
    line( ytim, yhat( :, maxidx ), 'color', colors( 3, : ) );
    xpos                        = min( xlim ) + 0.9 * diff( xlim );
    ypos                        = min( ylim ) + 0.9 * diff( ylim );
    th                          = text( xpos, ypos ...
        , sprintf( 'tr%d; cc=%0.3g', maxidx, maxval ) ); 
    set( th, 'color', colors( 3, : ) )
    
    % input - analog
    subplot( ax( 2 ) );
    line( ytim, y1( :, maxidx ), 'color', colors( 2, : ) );
    axis tight
    local_calibration( [ 100 1 ], { 'ms', zstr } );
    line( xlim, [ 0 0 ], 'color', blackColor, 'linestyle', '--' )
  
    % filters
    subplot( ax( 3 ) )
    if ismember( ftype,  { 'cwiener', 'xwiener', 'csta', 'xsta', 'cwsta', 'xwsta' } )
        plot( ftim, f1, 'b', ftim, f1 + f1sem, '--b', ftim, f1 - f1sem, '--b' );
        title( sprintf( 'STA (%d spks; W=%0.3g ms)' ...
            , stats.nspikes, stats.widsta / stats.Fs / 2 ) )
    else
        plot( ftim, f, 'b', ftim, f + fsem, '--b', ftim, f - fsem, '--b' );
        title( sprintf( '%s (%d spks; W=%0.3g ms)' ...
            , upper( ftype ), stats.nspikes, stats.widsta / stats.Fs / 2 ) )
    end
    axis tight
    if ismember( ftype,  { 'cwiener', 'xwiener', 'csta', 'xsta', 'cwsta', 'xwsta' } )
        subplot( ax( 4 ) )
        imagesc( ftim, 1 : nibins, f' ), axis xy
        set( ax( 4 ), 'ytick', ( 1 : nibins ) - 0.5, 'yticklabel', bedges( 1 : nibins ) / Fs )
        hold on
        plot( -mean( bedges( [ 1 : end - 1; 2 : end ] ) ) / Fs * 2, ( 1 : nibins ), '.-k' ) % twice the mean ISI
        line( [ 0 0 ], ylim, 'color', whiteColor, 'linestyle', '--' )
        title( sprintf( '%s (%d-%d spks)', [ ftype( 1 ) upper( ftype( 2 : end ) ) ] ...
            , min( stats.cpb ), max( stats.cpb ) ) )
    else
        ax( 4 )                 = -ax( 4 );
    end
    
    % filter optimization (time lags)
    if length( lagsr ) == 1
        ax( [ 7 8 ] )           = -ax( [ 7 8 ] );
    else
        subplot( ax( 7 ) )
        imagesc( lags / Fs, ftim, cci ), axis xy
        title( sprintf( 'CC; optimal lag: %0.2g', optlag / Fs ) )
        line( optlag / Fs * [ 1 1 ], ylim, 'color', redColor, 'linestyle', '--' )
        
        subplot( ax( 8 ) )
        plot( lags / Fs, max( cci, [], 1 ) )
        title( 'max CC/lag' )
    end

    % input-reconstruction CC
    subplot( ax( 5 ) )
    line( ftc, cct, 'color', grayColor )
    axis tight
    ylims                       = ylim;
    sc                          = ( act - min( act ) ) ...
        / ( max( act ) - min( act ) ) * ( max( ylims ) - min( ylims ) ) + min( ylims );
    line( ftc, sc, 'color', purpleColor, 'linewidth', 2 )
    line( ftc, cct, 'color', grayColor )
    line( ftc, cct( :, maxidx ), 'color', colors( 3, : ) );
    line( ftc, nanmean( cct, 2 ), 'color', blackColor, 'linewidth', 2 );
    xlim( [ min( ftc ) max( ftc ) ] )
    title( sprintf( '%s: %0.2g (p=%0.2g)', ccstr, ccAll, stats.pval ) )
    
    % scattergram
    subplot( ax( 6 ) )
    plot( ysct( :, 1 ), ysct( :, 2 ), '.k' )
    axis tight
    title( sprintf( '%0.3g', ccLinear ) )
    xlabel( 'yhat', 'color', colors( 3, : ) )
    ylabel( '<y>|yhat', 'color', colors( 2, : ) )
    set( ax( 6 ), 'yAxisLocation', 'right' )

    % coherence
    subplot( ax( 9 ) )
    imagesc( frq, 1 : nt, cohs' ), axis xy
    mcoh                        = nanmean( cohs, 2 );
    [ maxcoh, maxidx ]          = max( mcoh ); 
    ylims                       = ylim;
    sc                          = ( mcoh - min( mcoh ) ) ...
        / ( max( mcoh ) - min( mcoh ) ) * ( max( ylims ) - min( ylims ) ) + min( ylims );
    line( frq, sc, 'color', colors( 2, : ), 'linewidth', 2 );
    xlabel( 'Freq [Hz]' )
    title( sprintf( '<C> @ %0.2g Hz = %0.2g', frq( maxidx ), maxcoh ) )

    % organize
    ax( 10 )                    = -ax( 10 );
    for i                       = 1 : length( ax )
        if ax( i ) < 0
            subplot( -ax( i ) )
            axis off
            continue
        end
        subplot( ax( i ) );
        set( gca, 'tickdir', 'out', 'box', 'off' )
        if i == 1 || i == 2
            set( gca, 'ytick', [] )
            xlim( [ min( ytim ) max( ytim ) ] )
            axis off
        end
        if ismember( i, [ 3 5 6 7 ] )
            line( xlim, [ 0 0 ], 'color', blackColor, 'linestyle', '--' )
            line( [ 0 0 ], ylim, 'color', blackColor, 'linestyle', '--' )
        elseif i == 8
            line( [ 0 0 ], ylim, 'color', blackColor, 'linestyle', '--' )
        end
        
    end
    local_fig_title( tstr );
    colormap( myjet )
end

return % streconstruct

%------------------------------------------------------------------------
% h = local_shift( p )
% columns of a matrix by integer numbers
%------------------------------------------------------------------------
function [ y, idx0, idx1 ]      = local_shift( x, nshift, padwith )

nargs                           = nargin;
if nargs < 1 || isempty( x )
    return
end
if ~ismatrix( x )
    return
end
if nargs < 2 || isempty( nshift )
    nshift                      = 0;
end
if any( nshift ~= round( nshift ) )
    return
end
[ m, n ]                        = size( x );
if length( nshift ) == 1
    nshift                      = nshift * ones( 1, n );
elseif length( nshift ) ~= n
    return
end
if nargs < 3 || isempty( padwith )
    padwith                     = NaN;
end

idx                             = ( 1 : m )' * ones( 1, n );
idx                             = idx - ones( m, 1 ) * nshift( : ).';
idx( idx <= 0 | idx > m )       = 0;
idxhat                          = idx + m * ones( m, 1 ) * ( 0 : ( n - 1 ) );
idx1                            = find( idx > 0 );
idx0                            = idxhat( idx1 );
y                               = padwith * ones( m, n );
y( idx1 )                       = x( idx0 );

return % local_shift

%------------------------------------------------------------------------
% y = local_rankcols( x )
% rank matrix column-wise
%------------------------------------------------------------------------
function y = local_rankcols( x )

[ m, n ]                        = size( x );
if m == 1
    x                           = x';
    m                           = n;
    n                           = 1;
end
nans                            = isnan( x );
ridx                            = m : -1 : 1;
cidx                            = 1 : n;
[ ~, idx ]                      = sort( [ x x( ridx, : ) ] );
[ ~, ranks ]                    = sort( idx );
y                               = ( ranks( :, cidx ) + ranks( ridx, cidx + n ) ) / 2;
y( nans )                       = NaN;

return % local_rankcols

%------------------------------------------------------------------------
% Y = local_firfilt( x, W )
% zero-phase lag low-pass filtering of x's columns with the FIR W
%------------------------------------------------------------------------
function Y = local_firfilt( x, W )

C                               = length( W );
D                               = ceil( C / 2 ) - 1;
Y                               = filter( W, 1, [ flipud( x( 1 : C, : ) ); x; flipud( x( end - C + 1 : end, : ) ) ] );
Y                               = Y( 1 + C + D : end - C + D, : );

return % local_firfilt

%------------------------------------------------------------------------
% idx = local_make_ai( n )
% make auto-indices for vector of n elements
%------------------------------------------------------------------------
function idx = local_make_ai( n )

m                               = n * ( n - 1 ) / 2;
i                               = zeros( m, 1 );
j                               = ones( m, 1 );
p                               = ( n - 1 ) : -1 : 2;
i( cumsum( [ 1 p ] ) )          = 1;
i                               = cumsum( i );
j( cumsum( p ) + 1 )            = 2 - p;
if n > 1
    j( 1 )                      = 2;
end
j                               = cumsum( j );
idx                             = [ i j ];

return % local_make_ai

%------------------------------------------------------------------------
% blocks = local_makeblocks( n, blocksize )
% for batch processing
%------------------------------------------------------------------------
function blocks = local_makeblocks( n, blocksize )

nblocks                         = ceil( n / blocksize );
blocks                          = [ 1 : blocksize : blocksize * nblocks; blocksize : blocksize : blocksize * nblocks ]';
blocks( nblocks, 2 )            = n;
blocks( blocks > n )            = n;

return % local_makeblocks

%------------------------------------------------------------------------
% [ sval, eval, idx ] = local_make_equal_bins( x, n )
% assign indices to data s.t. equally-populated bins
%------------------------------------------------------------------------
function [ sval, eval, idx ]    = local_make_equal_bins( x, n )

% prepare
x                               = x( : );
sx                              = size( x );
nnans                           = ~isnan( x );
lx                              = sum( nnans );

% determine edges
[ xs, ix ]                      = sort( x );
bsize                           = lx / n;
si                              = cumsum( [ ones( 1, sx( 2 ) ); ones( n - 1, 1 ) * bsize ], 1 );
ei                              = cumsum( ones( n, 1 ) * bsize, 1 );
ios                             = ones( n, 1 ) * ( 0 : sx( 1 ) : ( sx( 1 ) * ( sx( 2 ) - 1 ) ) );
si                              = floor( si + ios );
ei                              = ceil( ei + ios );
si( si < 1 )                    = 1;
ei( ei > lx )                   = lx;
sval                            = xs( si );
eval                            = xs( ei );

% identicals
rmv                             = find( diff( sval ) == 0 );
if ~isempty( rmv )
    sval( rmv )                 = [];
    eval( rmv )                 = [];
end
rmv                             = find( diff( eval ) == 0 );
if ~isempty( rmv )
    sval( rmv )                 = [];
    eval( rmv )                 = [];
end

% actually bin
idx                             = zeros( sx );
for i                           = 1 : n
    for c                       = 1 : sx( 2 )
        idx( ix( si( i, c ) : ei( i, c ) ), c ) = i;
    end
end

return % local_make_equal_bins

%------------------------------------------------------------------------
% [ y, f ] = local_mtcsd1( x, nFFT, Fs, WinLength, nOverlap, NW, Detrend )
% cross-spectra between one signal (first column) and all others
%------------------------------------------------------------------------
function [ y, f ] = local_mtcsd1( x, nFFT, Fs, WinLength, nOverlap, NW, Detrend )

nTapers                         = 2 * NW - 1;
winstep                         = WinLength - nOverlap;
nChannels                       = size( x, 2 );
nSamples                        = size( x, 1 );

% check for column vector input
if nSamples == 1
	x                           = x';
	nSamples                    = size( x, 1 );
	nChannels                   = 1;
end

% calculate number of FFTChunks per channel
nFFTChunks                      = 1 + floor( ( ( nSamples - WinLength ) / winstep ) );

% allocate memory
y                               = complex( zeros( nFFT, 2, nChannels ) );   % output array
Periodogram                     = complex( zeros( nFFT, nTapers, nChannels ) ); % intermediate FFTs
Temp1                           = complex( zeros( nFFT, nTapers ) );
Temp2                           = complex( zeros( nFFT, nTapers ) );
Temp3                           = complex( zeros( nFFT, nTapers ) );
eJ                              = complex( zeros( nFFT, 1 ) );

% calculate Slepian sequences
Tapers                          = dpss( WinLength, NW, nTapers, 'calc' );

% compute tapered periodogram with FFT 
TaperingArray                   = repmat( Tapers, [ 1 1 nChannels ] );
for j                           = 1 : nFFTChunks
	Segment                     = x( ( j - 1 ) * winstep + ( 1 : WinLength ), : );
	if ~isempty( Detrend )
		Segment                 = detrend( Segment, Detrend );
    end
	SegmentsArray               = permute(repmat(Segment, [ 1 1 nTapers ] ), [ 1 3 2 ] );
	TaperedSegments             = TaperingArray .* SegmentsArray;

	Periodogram( :, :, : )      = fft( TaperedSegments, nFFT );

	% products 
    Temp1                       = squeeze( Periodogram( :, :, 1 ) );        % Ch1
    for Ch2                     = 1 : nChannels 
        
        Temp2                   = squeeze( Periodogram( :, :, Ch2 ) );
        Temp2                   = conj( Temp2 );
        Temp3                   = Temp1 .* Temp2;
        eJ                      = sum( Temp3, 2 );
        y( :, 2, Ch2 )          = y( :, 2, Ch2 ) + eJ / ( nTapers * nFFTChunks ); % cross

        Temp3                   = squeeze( Periodogram( :, :, Ch2 ) ) .* Temp2;
        eJ                      = sum( Temp3, 2 );
        y( :, 1, Ch2 )          = y( :, 1, Ch2 ) + eJ / ( nTapers * nFFTChunks ); % auto
    end
    
end

% select subset of y, set up f array
if ~any( any( imag( x ) ) )    % x purely real
	if rem( nFFT, 2 )
		select                  = 1 : ( nFFT + 1 ) / 2;
	else
		select                  = 1 : nFFT / 2 + 1;
	end
	y                           = y( select, :, : );
else
	select                      = 1 : nFFT;
end
f                               = ( select - 1 )' * Fs / nFFT;

return % local_mtcsd1

%------------------------------------------------------------------------
% local_calibration( len, s )
% add x-y calibration bars
%------------------------------------------------------------------------
function local_calibration( len, s )

loc0                            = [ 0.1 0.1 ]; 
xlims                           = xlim;
ylims                           = ylim;
dx                              = diff( xlims );
dy                              = diff( ylims );
xy                              = [ min( xlims ) + loc0( 1 ) * dx min( ylims ) + loc0( 2 ) * dy ];
dxy0                            = [ dx dy ];
dxy( len >= 0 )                 = len( len >= 0 );
dxy( len < 0 & len >= -1 )      = dxy0( len < 0 & len >= -1 ) .* abs( len( len < 0 & len >= -1 ) );
dxy( len < -1 )                 = NaN;

x                               = xy( 1 ) + [ 0 0 dxy( 1 ) ];
y                               = xy( 2 ) + [ dxy( 2 ) 0 0 ];
line( x, y, 'color', [ 0 0 0 ], 'linewidth', 2 );
str1                            = sprintf( '%0.3g %s', dxy( 1 ), s{ 1 } );
str2                            = sprintf( '%0.3g %s', dxy( 2 ), s{ 2 } );
th( 1 )                         = text( xy( 1 ) + dxy( 1 ) / 2, xy( 2 ) - dxy( 2 ) / 2, str1 );
set( th( 1 ), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'rotation', 0 )
th( 2 )                         = text( xy( 1 ) - dxy( 1 ) / 2, xy( 2 ) + dxy( 2 ) / 2, str2 );
set( th( 2 ), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'rotation', 90 )

return % local_calibration

%------------------------------------------------------------------------
% th = local_fig_title( tstr )
% title for the figure
%------------------------------------------------------------------------
function local_fig_title( tstr )

idx                             = strfind( tstr, '_' ); 
idx                             = [ idx length( tstr ) + 1 ];
str                             = tstr( 1 : idx( 1 ) - 1 );
for i                           = 1 : ( length( idx ) - 1 )
    str                         = sprintf( '%s%s%s', str, '\_' ...
        , tstr( ( idx( i ) + 1 ) : ( idx( i + 1 ) - 1 ) ) ); 
end

ah0                             = gca;
axes( 'position', [ 0.5 0.95 0.01 0.01 ] );
axis off;
th                              = text( 0.5, 0.5, str, 'Fontsize', 12 );
set( th, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' )
subplot( ah0 );

return % local_fig_title

% EOF
