% stAnalysis            computational routine for the heart of the gwnAnalysis
%
% call                  [ stats, details, yhat, fig ] = stAnalysis( st, wn )
%
% gets                  st              sparse matrix, trials in columns
%                       wn              full matrix, same size as st
%
% optional arguments (given as name/value pairs) passed to streconstruct:
%
%                       stHat           {[]}                    argument passed as x2 to streconstruct
%                       nfiltSEC        {0.05}      [s]         use a filter/bank of this half-width
%                       ftype           {'wiener'}              other options include 'sta', 'corr'; 'csta', 'cwiener'; 'xsta', 'xwiener'
%                       maxISIbins      {10}                    maximal number of ISI bins (relevant only for 2nd order models, 'csta', 'cwiener', 'xsta', 'xwiener')
%                       Fs              {2500}      [Hz]        frequency at which st and wn are sampled
%                                                               used to convert nfiltSEC to samples before passing to streconstruct
%                       prepMode        {'flat'}                used (locally and) in streconstruct
%                       minSpikesPerTrial {1}                   the minimal mean number of spikes/trial
%                       minTrials       {4}                     the minimal number of trials
%                       nreps0          {10}                    shuffles (for reconstruction quality)
%                       doRel           {0}                     argument passed to streconstruct
%
% optional arguments (given as name/value pairs) passed to streconstructXval:
%
%                       cvMode          {10}                    K-fold cross-validation. 1 indicates jackknife;
%                                                                 other numbers indicate the number of folds
%
% optional arguments (given as name/value pairs) used locally for computing R, R-profile, and P3:
%
%                       doRelLocal      {-1}                    1, compute R-profile and P3 (if Q is sign.)
%                                                               -1, compute only R (anyhow)
%                       pTH             {0.05}                  alpha level; intersected with doRelLocal (R profile, P3)
%                                                                       and doTiming (P, P2)
%                       tt              {see@left}  [samples]   half-window sizes used for convolution (R-profile, P3) and jittering (precision: P, P2)
%                                                               [ 1 2 4 8 12 16 24 32 64 128 256 ]
%                       ftypeR          {'gauss'}               'ftype' argument of streliability
%                       doFocusedR      {1}                     'doFocused' argument of streliability
%
% optional arguments (given as name/value pairs) used locally for computing P, P2:
%
%                       doTiming        {1}         logical     compute precision (using stmix and stahat) locally (if Q is sign.)
%                                                               using 2 methods, both based on spike jittering: P: fixed filter; P2: data-dependent filter
%                       nreps1          {0}                     jitter iterations (for precision by fixed filter deterioration)
%                       nreps2          {0}                     jitter iterations (for precision by data-dependent filter deterioration)
%                       jmode           {1}                     jitter mode, passed to stmix; 1: spike-time jitter; 2: interval jitter
%
% optional arguments (given as name/value pairs) used locally (flow control)
%:
%                       toIgnoreSEC     {0.1}       [s]         used locally for input selection (e.g. ignore the first part)
%                       verboseCV       {0}         flag        during cross-validation (see step 1 below)
%                       verboseTim      {0}         flag        during timing computations (see step 4 below)
%                       graphics        {0}                     1: plots; -1: plots only if sig. Q
%
% returns               stats           summary structure with fields:
%                                           Q       quality vector: [ mean SEM pval ]
%                                           R       reliability vector: [ mean SEM pval ]
%                                           P       precision by fixed filter [ms]
%                                           T       time lag, [ms]
%                                           Rprof   [ mean SEM ] at each tt
%                                           Rp      [ SD mean SEM ] value of Rprof at peak (may be different from all tt's; SD in ms)
%                                           P2      precision by data-dependent filter, [ms]
%                                           MI      [bits]; MI details: [ debiased MI; MI (biased); bias; stimulus entropy; reconsuction entropy; joint entropy (biased) ]
%                                           Coh     coherence in the bandpass (0-40 Hz)
%
%                       details         detailed summary structure
%                                       based second output of streconstruct, and augmented with
%                                           tt      actual half-windows used as SDs and for jittering, [samples]
%
%                                           ccAll   mean and SEM of cross-validated Q
%                                           pval    pval of cross-validated Q
%                                           cc      trial-specific Q
%                                           pvals   trial-specific p-value
%
%                                           widAll  FWHH for mean cross-correlations between cross-validated yhat and the input
%
%                                   Reliability
%                                           ttrelF  list of SDs for focused computation of reliability profile
%                                           RprofF  Rprof for the focused SDs (ttrelF)
%
%                                   Precision (fixed filter)
%                                           ccIn    rho for each SD in tt: [ mean SEM ]
%                                           pIn     p-value for each rho in tt
%                                           ttInF   list of jittering window sizes for focused computation of precision
%                                           ccInF   same as ccIn, for focused search
%                                           pInF    same as pIn, for focused search
%
%                                   Precision (data dependent filter)
%                                           ccEx    same as ccIn, for data-dependent filter
%                                           pEx     same as pIn, for data-dependent filter
%                                           ttExF   same as ttInF, for data-dependent filter
%                                           ccExF   same as ccInF, for data-dependent filter
%                                           pExF    same as pInF, for data-dependent filter
%                                           pExJit
%                                           pExJitF
%
%                       yhat            wn reconstruction, based on streconstruct
%                       fig             handle to figure
%
% does
% (1) --optional--  computes cross-validated filters, uses that for Q-based stats       cvMode      streconstructXval.m
% (2)               computes filters, uses that for T stats                                         streconstruct.m
% (3) --optional--  computes reliability and R-profile                                  doRelLocal  streliability
% (4) --optional--  computes precision (P, P2) by spike jittering                       doTiming    stmix, stahat
% (5) --optional--  summarizes graphically                                              graphics
% (6) --optional--  computes MI and coherence                                           mutinf, local_mtcsd1
%
% calls:                ParseArgPairs
%                       calc_bias_PT, calc_fwhh, calc_spearman
%                       isisort, mutinf, my_xcorr, nangeomean
%                       stahat, stmix, streliability, streconstruct, streconstructXval, stcodertypes
%                       zsc (mex)

% 07-jan-15 ES

% last update
% August 2022

function [ stats, details, yhat, fig ] = stAnalysis( st, wn, varargin )

%----------------------------------------------------------------------
% constants
%----------------------------------------------------------------------
tt_DFLT                         = [ 1 2 4 8 12 16 24 32 64 128 256 ];       % [samples] - jitter time scales (half windows)
minQ                            = 0;                                        % or 0.02

% MI constants
MI_nbins                        = 20;                                       % argument to bindata
MI_debias                       = 'bayes';                                  % see calc_bias_PT

% Coherence constants
Coh_fROI                        = [ 0 200 ];
Coh_M                           = 0.25;
Coh_mtNW                        = 3;
Coh_dflag                       = '';

% graphics
colors                          = [ 0 0 0.7; 0 0 0; 1 0 0 ];                % x, y, yhat
blackColor                      = [ 0 0 0 ];
grayColor                       = [ 1 1 1 ] * 0.6;
purpleColor                     = [ 1 0.5 1 ];
xscale                          = 'log';
pTHplot                         = 0.05;                                     % Qpval cross val
USF                             = 10;                                       % argument to imupsample for joint PDF

%----------------------------------------------------------------------
% defaults
%----------------------------------------------------------------------
Coh_fROI_stats_DFLT             = [ 0 40 ];                                 % bandpass to compute mean coherence at

%----------------------------------------------------------------------
% initialize
%----------------------------------------------------------------------
stats                           = [];
details                         = [];
yhat                            = [];
fig                             = [];

R                               = NaN;
Rprof                           = NaN;
P                               = NaN;
P2                              = NaN;
Rp                              = NaN;

%----------------------------------------------------------------------
% arguments
%----------------------------------------------------------------------
nargs                           = nargin;
if nargs < 2 || isempty( st ) || isempty( wn )
    return
end
if size( st, 1 ) ~= size( wn, 1 )
    fprintf( '%s: input size mismatch!', upper( mfilename ) )
    return
end
nt                              = size( st, 2 );
if nt > 1 && size( wn, 2 ) == 1
    wn                          = wn * ones( 1, nt );
end
if size( wn, 2 ) ~= nt
    fprintf( '%s: input size mismatch!', upper( mfilename ) )
    return
end

[ doRelLocal, stHat ...
    , Fs ...
    , toIgnoreSEC, nfiltSEC, ftype, wfvar, maxISIbins ...
    , prepMode, pwmode ...
    , minTrials, minSpikesPerTrial ...
    , doTiming, doRel, pTH, cvMode ...
    , ftypeR, doFocusedR ...
    , nreps0, nreps1, nreps2, jmode ...
    , tt ...
    , doMI, doCoh, Coh_fROI_stats ...
    , verboseCV, verboseTim ...
    , graphics ]                = ParseArgPairs(...
    { 'doRelLocal', 'stHat' ...
    , 'Fs' ...
    , 'toIgnoreSEC', 'nfiltSEC', 'ftype', 'wfvar', 'maxISIbins' ...
    , 'prepMode', 'pwmode' ...
    , 'minTrials', 'minSpikesPerTrial' ...
    , 'doTiming', 'doRel', 'pTH', 'cvMode' ...
    , 'ftypeR', 'doFocusedR' ...
    , 'nreps0', 'nreps1', 'nreps2', 'jmode' ...
    , 'tt' ...
    , 'doMI', 'doCoh', 'Coh_fROI_stats' ...
    , 'verboseCV', 'verboseTim' ...
    , 'graphics' }...
    , { -1, [] ...
    , 2500 ...
    , 0.1, 0.05, 'wiener', 0.9, 10 ...
    , 'flat', 'pca' ...
    , 4, 1 ...
    , 1, 0, 0.05, 10 ...
    , 'gauss', 1 ...
    , 10, 10, 0, 1 ...
    , tt_DFLT ...
    , 1, 1, Coh_fROI_stats_DFLT ...
    , 0, 0 ...
    , [ 0 1 ] }...
    , varargin{ : } );
cvMode                          = cvMode( 1 );
if cvMode < 0 || cvMode ~= round( cvMode )
    cvMode                      = 0;
end
nfilt                           = round( nfiltSEC * Fs );
maxlag                          = max( abs( nfilt ) );
switch prepMode
    case { 'flat', 'unbiased' }
        idx0                    = floor( max( toIgnoreSEC - max( abs( nfiltSEC ) ), 0 ) * Fs ) + 1;
    case 'pad'
        idx0                    = floor( toIgnoreSEC * Fs ) + 1;
    case 'none'
        idx0                    = 1;
    otherwise
        error( 'unsupported prepMode' )
end
if isequal( ftype( 1 : 2 ), 'pw' )
    prew                        = 1;
    ftype                       = ftype( 3 : end );
else
    prew                        = 0;
end

%----------------------------------------------------------------------
% preps
%----------------------------------------------------------------------
% ignore beginning of data:
idx                             = idx0 : size( wn, 1 );
nsamples                        = length( idx );
if isempty( stHat )
    stHat                       = sparse( size( st, 1 ), size( st, 2 ) );
end
% prewhiten:
if prew
    same_wn                     = isequal( wn( :, 1 ) * ones( 1, nt - 1 ), wn( :, 2 : nt ) );
    if same_wn
        z                       = whitensig( wn( idx, 1 ), 'wmode', pwmode, 'wfvar', wfvar );
        wn( idx, : )            = z * ones( 1, nt );
    else
        z                       = whitensig( wn( idx, : ), 'wmode', pwmode, 'wfvar', wfvar );
        wn( idx, : )            = z;
    end
end
% preps for graphics
if any( graphics ) && doTiming
    xlims                       = [ min( tt ) max( tt ) ] / Fs * 1000;
    xticks                      = tt / Fs * 1000;
end

%----------------------------------------------------------------------
% step 1. model (+cross-validation):
%----------------------------------------------------------------------
if cvMode
    [ yhatCV, ccCV, pvalCV, models ] = streconstructXval( st( idx, : ), wn( idx, : ), stHat( idx, : ) ...
        , 'nfolds', cvMode, 'nrepsCV', nreps0 ...
        , 'nfilt', nfilt, 'ftype', ftype, 'maxISIbins', maxISIbins, 'Fs', Fs / 1000, 'lagsr', 0 ...
        , 'prepMode', prepMode ...
        , 'minSpikesPerTrial', minSpikesPerTrial, 'minTrials', minTrials ...
        , 'nreps', 0, 'doRel', 0, 'graphics', 0, 'verbose', verboseCV );
    [ ~, midx ]                 = sort( [ models.testset ] );
    ccmix                       = [ models.cc ];
    pvalsmix                    = [ models.pval ];
    ccsCV                       = ccmix( midx );
    pvalsCV                     = pvalsmix( midx );
    if isempty( ccsCV )
        ccsCV                   = NaN * ones( 1, nt );
    end
    if isempty( pvalsCV )
        pvalsCV                 = NaN * ones( 1, nt );
    end
end

% summary (f=<f(CV)>)
if graphics( 1 ) == 1 || ( graphics( 1 ) == -1 && pvalCV <= pTHplot )
    graphics1                   = 1;
else
    graphics1                   = 0;
end
[ yhat, stats, ~, ax ]          = streconstruct( st( idx, : ), wn( idx, : ), 'x2', stHat( idx, : ) ...
    , 'nfilt', nfilt, 'ftype', ftype, 'maxISIbins', maxISIbins, 'Fs', Fs / 1000, 'lagsr', 0 ...
    , 'prepMode', prepMode ...
    , 'minSpikesPerTrial', minSpikesPerTrial, 'minTrials', minTrials ...
    , 'nreps', nreps0 * double( cvMode == 0 ), 'doRel', doRel, 'graphics', graphics1 );

if isempty( stats )
    return
end
stats.tt                        = tt;

% fix by the results of the cross-validated runs
if cvMode
    stats.ccAll                 = ccCV;
    stats.pval                  = pvalCV;
    yhat                        = yhatCV;
    stats.cc                    = ccsCV;
    stats.pvals                 = pvalsCV;
    rwn                         = local_rankcols( wn( idx, : ) );
    cct                         = my_xcorr( local_rankcols( yhat ), rwn, maxlag, -1 );
    stats.widAll                = calc_fwhh( nanmean( cct, 2 ) );
end

% initial graphics
if graphics1
    fig                         = gcf;
    if cvMode
        ftim                    = stats.ftim;
        act                     = mean( my_xcorr( rwn, rwn, maxlag, -1 ), 2 );
        [ maxval, maxidx ]      = max( stats.cc );
        ytim                    = ( 1 : nsamples )' / Fs * 1000;
        
        subplot( ax( 2 ) )
        cla
        line( ytim, zsc( yhat( :, maxidx ), 1 ), 'color', colors( 3, : ) );
        xpos                    = min( xlim ) + 0.9 * diff( xlim );
        ypos                    = min( ylim ) + 0.9 * diff( ylim );
        th                      = text( xpos, ypos, sprintf( 'tr%d; q=%0.2g', maxidx, maxval ) );
        set( th, 'color', colors( 3, : ) )
        line( ytim, zsc( wn( idx, maxidx ), 1 ), 'color', colors( 2, : ) );
        axis tight
        local_calibration( [ 100 1 ], { 'ms', 'Z' } );
        line( xlim, [ 0 0 ], 'color', blackColor, 'linestyle', '--' )
        
        axis off
        
        subplot( ax( 5 ) )
        cla
        title( sprintf( 'Q = %0.2g (p=%0.2g)', stats.ccAll( 1 ), stats.pval ) )
        line( ftim, cct, 'color', grayColor )
        axis tight
        ylims                   = ylim;
        sc                      = ( act - min( act ) ) / ( max( act ) - min( act ) ) * ( max( ylims ) - min( ylims ) ) + min( ylims );
        line( ftim, sc, 'color', purpleColor, 'linewidth', 2 )
        line( ftim, cct, 'color', grayColor )
        line( ftim, cct( :, maxidx ), 'color', colors( 3, : ) );
        line( ftim, nanmean( cct, 2 ), 'color', blackColor, 'linewidth', 2 );
        xlim( [ min( ftim ) max( ftim ) ] )
        %alines( 0, 'y', 'color', blackColor, 'linestyle', '--' );
        %alines( 0, 'x', 'color', blackColor, 'linestyle', '--' );
        line( xlim, [ 0 0 ], 'color', blackColor, 'linestyle', '--' )
        line( [ 0 0 ], ylim, 'color', blackColor, 'linestyle', '--' )
        
    end
end

% summarize:
T                               = stats.lagsta / Fs * 1000;                 % [ms]; time lag, based on the STA filter
Q                               = stats.ccAll;
Qpval                           = stats.pval;

% initialize additional fields
ccIn                            = NaN * ones( length( tt ), 2 );
pIn                             = NaN * ones( length( tt ), 1 );
ccEx                            = ccIn;
pEx                             = pIn;
pExJit                          = pIn;

ttrelF                          = [];
RprofF                          = [];

ttInF                           = [];
ccInF                           = [];
pInF                            = [];

ttExF                           = [];
ccExF                           = ccInF;
pExF                            = pInF;
pExJitF                         = pInF;

MI                              = NaN * ones( 1, 6 );
mcoh                            = NaN;

nFFT                            = 2^floor( log2( Coh_M * Fs ) );
frq                             = ( 0 : nFFT / 2 )' * Fs / nFFT;            % all frequencies up to Nyquist
frqidx                          = frq >= min( Coh_fROI ) & frq <= max( Coh_fROI );
frq                             = frq( frqidx );
cohs                            = NaN( sum( frqidx ), size( yhat, 2 ) );
phs                             = cohs;

%----------------------------------------------------------------------
% step 2. reliability (+precision method 3; based on cross-trial correlation coefficient)
%----------------------------------------------------------------------
if Qpval <= pTH && Q( 1 ) >= minQ && doRelLocal || doRelLocal == -1         % do this anyhow.. doesn't waste much time if doRelLocal is -1
    
    if graphics1
        subplot( abs( ax( 10 ) ) );
        cla
    end
    if doRelLocal == 1
        [ R, vrreljit, vrreljitF, ttrelF ] = streliability( st( idx, : ), yhat, 's', tt ...
            , 'ftype', ftypeR, 'doFocused', doFocusedR, 'graphics', graphics( 1 ) );
        Rprof                   = [ mean( vrreljit, 2 ) std( vrreljit, [], 2 ) / sqrt( size( vrreljit, 2 ) ) ];
        RprofF                  = [ mean( vrreljitF, 2 ) std( vrreljitF, [], 2 ) / sqrt( size( vrreljitF, 2 ) ) ];
        
        if isnan( RprofF )
            RprofF               = [];
        end
        
        % compute peak of all tested SDs
        RprofU                  = [ Rprof; RprofF ];
        ttU                     = [ tt( : ); ttrelF ];
        [ ~, sidx ]             = unique( ttU );
        mat                     = [ ttU( sidx ) RprofU( sidx, : ) ];
        [ ~, maxidx ]           = max( mat( :, 2 ) );
        Rp                      = mat( maxidx, : );
        Rp( 1 )                 = Rp( 1 ) / Fs * 1000;                      % [samples] -> [ms]
        
        % call stcodertypes
        mm                      = mat( :, 2 );
        ss                      = mat( :, 3 );
        nn                      = nt;
        bb                      = mat( :, 1 )/ Fs * 1000;
        rFull0                  = [ vrreljit; vrreljitF ];
        rFull                   = rFull0( sidx, : );
        codeType                = stcodertypes( mm, ss, nn, bb, rFull );
        
    elseif doRelLocal == -1
        
        R                       = streliability( st( idx, : ), yhat, 's', NaN ...
            , 'ftype', 'gauss', 'graphics', 0 );
        
    end
    
    if graphics1
        xticks3l                = get( gca, 'xticklabel' );
        xticks3                 = NaN( 1, size( xticks3l, 1 ) );
        for tn                  = 1 : size( xticks3l, 1 )
            xticks3( tn )       = str2double( xticks3l( tn, : ) );
        end
        set( gca, 'xticklabel', xticks3 / Fs * 1000 )
        title( sprintf( 'Rel = %0.2g (p=%0.2g)', R( 1 ), R( 3 ) ) )
        xlabel( '\sigma [ms]' )
    end
    
else
    
    codeType                    = NaN;
    
end

%----------------------------------------------------------------------
% step 3. temporal precision:
%----------------------------------------------------------------------
if doTiming && Qpval <= pTH && Q( 1 ) >= minQ
    
    % method 1 - internal (see how much the spikes need to be jittered to
    % kill the reconstruction with a FIXED filter)
    if nreps1 > 0
        if verboseTim
            fprintf( '\tPrecision (method1 - fixed model); jittering scale [samples]: ' )
        end
        
        nrepsK                  = nreps1;
        for k                   = 1 : length( tt )
            ttK                 = tt( k );
            if verboseTim
                fprintf( 1, '%d ', ttK )
            end
            st0                 = st( idx, : );
            st0( st0 > 1 )      = 1;
            xm                  = stmix( st0, ttK, nrepsK, jmode );         % for now, decimate
            xidx                = ones( nrepsK, 1 ) * ( 1 : nt ) ;
            xidx                = xidx( : );
            yhatJ               = NaN * ones( size( xm ) );
            for h               = 1 : length( models )
                if isempty( models( h ).testset )
                    continue
                end
                jix             = ismember( xidx, models( h ).testset );
                xt              = isisort( xm( :, jix ), models( h ).isiedges );
                yhatJ( :, jix ) = stahat( xt, models( h ).f, nsamples, models( h ).optlag );
            end
            ccJitK              = calc_spearman( yhatJ, wn( idx, xidx ) );
            ccIn( k, : )        = [ nanmean( ccJitK ) ...
                std( ccJitK, [], 2 ) / sqrt( size( ccJitK, 2 ) ) ];
            pIn(  k, : )        = ranksum( stats.cc( : ), ccJitK( : ) ...
                , 'tail', 'right' );                                        % is the reconstruction - at this jitter - worse than for raw data?
        end
        hwnum                   = find( pIn <= pTH, 1, 'first' );
        if ~isempty( hwnum )
            if tt( hwnum ) == 1
                P               = tt( hwnum ) / Fs * 1000;                  % [ms]
            else
                % repeat everything for the range between tt( hwnum ) and the one before it
                ttInF           = tt( hwnum - 1 ) : tt( hwnum );
                nFocused        = length( ttInF );
                ccInF           = NaN * ones( nFocused, 2 );
                pInF            = NaN * ones( nFocused, 1 );
                for k           = 1 : length( ttInF )
                    ttK         = ttInF( k );
                    st0         = st( idx, : );
                    st0( st0 > 1 )      = 1;
                    xm          = stmix( st0, ttK, nrepsK, jmode );
                    xidx        = ones( nrepsK, 1 ) * ( 1 : nt ) ;
                    xidx        = xidx( : );
                    yhatJ       = NaN * ones( size( xm ) );
                    for h       = 1 : length( models )
                        if isempty( models( h ).testset )
                            continue
                        end
                        jix     = ismember( xidx, models( h ).testset );
                        xt      = isisort( xm( :, jix ), models( h ).isiedges );
                        yhatJ( :, jix ) = stahat( xt, models( h ).f, nsamples, models( h ).optlag );
                    end
                    ccJitK      = calc_spearman( yhatJ, wn( idx, xidx ) );
                    ccInF( k, : )       = [ nanmean( ccJitK ) ...
                        std( ccJitK, [], 2 ) / sqrt( size( ccJitK, 2 ) ) ];
                    pInF(  k, : )       = ranksum( stats.cc( : ) ...
                        , ccJitK( : ), 'tail', 'right' );
                end
                hwnumF          = find( pInF <= pTH, 1, 'first' );          % could happen since nreps1 is not inf
                if isempty( hwnumF )
                    if hwnum > 1
                        hwnum 	= hwnum + [ -1 0 ];
                    end
                    P         	= nangeomean( tt( hwnum ) / Fs * 1000 );    % [ms]
                else
                    P           = ttInF( hwnumF ) / Fs * 1000;              % [ms]
                end
            end
        end
        
        if verboseTim
            fprintf( ': %0.3g ms\n', P  )
        end
        
    end
    
    
    % method 2 - external (see how much the spikes need to be jittered to
    % kill the correspondence with the input, i.e. build a NEW filter each time)
    if nreps2 > 0
        if verboseTim
            fprintf( 1, '\tPrecision (method2: jittering - data-dependent models):\n' )
        end
        for k                   = 1 : length( tt )
            if verboseTim
                fprintf( 1, '\t%d %s:\n', tt( k ), upper( mfilename ) )
            end
            ccExCV              = NaN * ones( nreps2, 2 );
            ppExCV              = NaN * ones( nreps2, 1 );
            ccJitK              = NaN * ones( nreps2, nt );
            for h               = 1 : nreps2
                sthat           = stmix( st( idx, : ), tt( k ), 1, jmode );
                if verboseTim
                    fprintf( 1, '\t\t' )
                end
                [ ~, ccExCV( h, : ), ppExCV( h, : ), modelsEx ] = streconstructXval( sthat, wn( idx, : ), [] ...
                    , 'nfolds', cvMode, 'nrepsCV', nreps0 ...
                    , 'minSpikesPerTrial', minSpikesPerTrial, 'minTrials', minTrials ...
                    , 'nfilt', nfilt, 'ftype', ftype, 'wfvar', wfvar, 'maxISIbins', maxISIbins ...
                    , 'Fs', Fs / 1000, 'lagsr', 0, 'rhalfwin', [], 'nreps' , 0 );
                ccjit           = [ modelsEx.cc ];
                ccJitK( h, 1 : length( ccjit ) ) = ccjit;
            end
            ccEx( k, : )        = nanmean( ccExCV );
            pExJit( k, : )      = nangeomean( ppExCV );                 	% is the x-val reconstruction - at this jitter - significant?
            pEx( k, : )         = ranksum( stats.cc( : ), ccJitK( : ), 'tail', 'right' );   % is the x-val reconstruction - at this jitter - worse than for raw data?
            if verboseTim
                fprintf( 1, '\n' )
            end
        end
        hwnum                   = find( pEx <= pTH, 1, 'first' );
        if ~isempty( hwnum )
            if tt( hwnum ) == 1
                P2              = tt( hwnum ) / Fs * 1000;                  % [ms]
            else
                % repeat everything for the range between tt( hwnum ) and the one before it
                ttExF           = tt( hwnum - 1 ) : tt( hwnum );
                nFocused        = length( ttExF );
                ccExF           = NaN * ones( nFocused, 2 );
                pExJitF         = NaN * ones( nFocused, 1 );
                pExF            = NaN * ones( nFocused, 1 );
                for k           = 1 : length( ttExF )
                    if verboseTim
                        fprintf( 1, '\t%d %s:\n', tt( k ), upper( mfilename ) )
                    end
                    ccExCV      = NaN * ones( nreps2, 2 );
                    ppExCV      = NaN * ones( nreps2, 1 );
                    ccJitK      = NaN * ones( nreps2, nt );
                    for h       = 1 : nreps2
                        sthat   = stmix( st( idx, : ), tt( k ), 1, jmode );
                        if verboseTim
                            fprintf( 1, '\t\t' )
                        end
                        [ ~, ccExCV( h, : ), ppExCV( h, : ), modelsEx ] = streconstructXval( sthat, wn( idx, : ), [] ...
                            , 'nfolds', cvMode, 'nrepsCV', nreps0 ...
                            , 'minSpikesPerTrial', minSpikesPerTrial, 'minTrials', minTrials ...
                            , 'nfilt', nfilt, 'ftype', ftype, 'wfvar', wfvar, 'maxISIbins', maxISIbins ...
                            , 'Fs', Fs / 1000, 'lagsr', 0, 'rhalfwin', [], 'nreps' , 0 );
                        ccjit   = [ modelsEx.cc ];
                        ccJitK( h, 1 : length( ccjit ) ) = ccjit;
                    end
                    ccExF( k, : )    = nanmean( ccExCV );
                    pExJitF( k, : )  = nangeomean( ppExCV );
                    pExF( k, : )     = ranksum( stats.cc( : ), ccJitK( : ), 'tail', 'right' );
                    if verboseTim
                        fprintf( 1, '\n' )
                    end
                end
                hwnumF          = find( pExF <= pTH, 1, 'first' );
                if isempty( hwnumF )                                        % could happen since nreps2 is not inf
                    if hwnum > 1
                        hwnum 	= hwnum + [ -1 0 ];
                    end
                    P2       	= nangeomean( tt( hwnum ) / Fs * 1000 );    % [ms]
                else
                    P2          = ttExF( hwnumF ) / Fs * 1000;              % [ms]
                end
            end
        end
        if verboseTim
            fprintf( '\tprecision = %0.3g ms\n', P2  )
        end
    end
    
    if graphics1
        cc0                     = Q;
        if nreps1 > 0
            ax( 7 )             = axes( 'position', [ 0.05  0.05 0.2 0.2 ] );
            subplot( ax( 7 ) )
            cla
            plot( tt / Fs * 1000, ccIn( :, 1 ), 'b', tt / Fs * 1000, ccIn( :, 1 ) + ccIn( :, 2 ), '--b' ...
                , tt / Fs * 1000, ccIn( :, 1 ) - ccIn( :, 2 ), '--b' );
            line( tt / Fs * 1000, ccIn( :, 1 ), 'Marker', '.', 'linewidth', 2, 'markerSize', 20 )
        end
        if nreps2 > 0
            ax( 8 )             = axes( 'position', [ 0.3   0.05 0.2 0.2 ] );
            subplot( ax( 8 ) )
            cla
            plot( tt / Fs * 1000, ccEx( :, 1 ), 'b', tt / Fs * 1000, ccEx( :, 1 ) + ccEx( :, 2 ), '--b' ...
                , tt / Fs * 1000, ccEx( :, 1 ) - ccEx( :, 2 ), '--b' );
            line( tt / Fs * 1000, ccEx( :, 1 ), 'Marker', '.', 'linewidth', 2, 'markerSize', 20 )
        end
        for aidx                = 7 : 8
            if aidx == 7 && nreps1 == 0 || aidx == 8 && nreps2 == 0
                continue
            end
            subplot( ax( aidx ) )
            set( gca, 'xscale', xscale, 'tickdir', 'out', 'box', 'off' )
            axis tight
            xlim( xlims )
            if isequal( xscale, 'log' )
                set( gca, 'xtick', xticks, 'xticklabel', xticks )
            end
            if aidx == 7
                title( sprintf( 'P = %0.2g ms', P ) )
                line( xlim, [ 0 0 ], 'color', blackColor, 'linestyle', '--' )
            else
                title( sprintf( 'P = %0.2g ms', P2 ) )
            end
            line( xlim, cc0( 1 ) * [ 1 1 ], 'color', blackColor, 'linestyle', '--' )
            line( xlim, ( cc0( 1 ) + cc0( 2 ) ) * [ 1 1 ], 'color', blackColor, 'linestyle', '--' )
            line( xlim, ( cc0( 1 ) - cc0( 2 ) ) * [ 1 1 ], 'color', blackColor, 'linestyle', '--' )
            ylims               = ylim;
            ylims( 2 )          = max( ylims( 2 ), ( cc0( 1 ) + cc0( 2 ) ) * 1.01 );
            ylim( ylims )
            xlabel( '\delta [ms]' )
            ylabel( 'CC' )
        end
    end % graphics
    
    if verboseTim
        fprintf( '\n\tDone!\n' )
    end
    
end % doTiming

%----------------------------------------------------------------------
% step 4. mutual information and coherence analyses
%----------------------------------------------------------------------
if doMI && Qpval <= pTH && Q( 1 ) >= minQ
    
    x                           = wn( idx, : );
    y                           = yhat;
    x                           = x( : );
    y                           = y( : );
    mm                          = [ min( x ) max( x ) ];
    bs                          = diff( mm ) / MI_nbins;
    xedges                      = ( mm( 1 ) - bs / 2 ) : bs : ( mm( 2 ) + bs / 2 );
    bdata                       = histcounts2( x, y, xedges, xedges );
    xbins                       = ( ( xedges( 1 : end - 1 ) + xedges( 2 : end ) ) / 2 )';
    ybins                       = xbins;
    px                          = sum( bdata, 2 );
    py                          = sum( bdata )';
    px                          = px / sum( px );
    py                          = py / sum( py );
    mi                          = mutinf( bdata );
    bias                        = calc_bias_PT( bdata, MI_debias );
    ixy                         = mi - bias;
    hx                          = local_calc_entropy( px );
    hy                          = local_calc_entropy( py );
    hxy                         = local_calc_entropy( bdata( : ) / sum( bdata( : ) ) );
    MI                          = [ ixy mi bias hx hy hxy ];
    if graphics1
        subplot( abs( ax( 8 ) ) )
        title( '' )
        cla
        axis off
        pos                     = get( gca, 'position' );
        ax( 8 )                 = axes( 'position', pos );
        imagesc( xbins, ybins, bdata' )
        axis xy
        xlabel( sprintf( 'y (%0.2g bits)', hx ) )
        ylabel( sprintf( 'yhat (%0.2g bits)', hy ) )
        set( gca, 'tickdir', 'out', 'box', 'off' )
        title( sprintf( 'MI = %0.2g bits', ixy ) )
        axis square
    end
    
end

if doCoh && Qpval <= pTH && Q( 1 ) >= minQ
    
    nWindow                     = nFFT;
    [ yo, frq ]                 = local_mtcsd1( [ wn( idx, 1 ) yhat ], nFFT, Fs, nWindow, nWindow/2, Coh_mtNW, Coh_dflag );
    yo                          = yo( frqidx, :, : );
    frq                         = frq( frqidx );
    p1                          = repmat( yo( :, 1, 1 ), [ 1 nt ] );
    p2                          = squeeze( yo( :, 1, 2 : ( nt + 1 ) ) );
    c12                         = squeeze( yo( :, 2, 2 : ( nt + 1 ) ) );
    cohs                        = single( abs( c12 .^ 2 ) ./ ( p1 .* p2 ) );
    phs                         = single( atan2( imag( c12 ), real( c12 ) ) );
    mcohs                       = mean( cohs, 2 );                          % mean over trials
    frqidxStats                 = frq >= min( Coh_fROI_stats ) & frq <= max( Coh_fROI_stats );
    mcoh                        = mean( mcohs( frqidxStats ) );             % mean over frequencies
    
    if graphics1 && ax( 4 ) < 0
        
        subplot( abs( ax( 4 ) ) )
        title( '' )
        cla
        axis off
        ax4                     = zeros( 2, 1 );
        ax4( 1 )                = axes( 'position', [ 0.3   0.3  0.2 0.1 ] );
        ax4( 2 )                = axes( 'position', [ 0.3   0.4  0.2 0.1 ] );
        
        subplot( ax4( 2 ) )
        plot( frq, mcohs, 'b' )
        set( gca, 'tickdir', 'out', 'box', 'off' )
        set( gca, 'xticklabel', '' )
        ylabel( 'Coherence' )
        tstr                    = sprintf( 'Coh = %2.2f', mcoh );
        ylims                   = ylim;
        xlims                   = xlim;
        th                      = text( diff( xlims ) * 0.75, diff( ylims ) * 0.75, tstr );
        set( th, 'HorizontalAlignment', 'Center' )
        
        subplot( ax4( 1 ) )
        uwphs                   = unwrap( local_circ_mean( phs' )' );
        if max( uwphs ) > 2 * pi
            uwphs               = uwphs - 2 * pi;
        end
        plot( frq, uwphs, 'b' )
        title( sprintf( '%0.2g rad', local_circ_mean( phs( frqidxStats ) ) ) )
        set( gca, 'tickdir', 'out', 'box', 'off' )
        line( xlim, [ 0 0 ], 'color', blackColor, 'linestyle', '--' )
        ylim( [ -1 1 ] * pi )
        xlabel( 'Frequency [Hz]' )
        ylabel( 'Phase [rad]' )
        
    end
    
end

%----------------------------------------------------------------------
% summarize in structure
%----------------------------------------------------------------------
stats.ccEx                      = ccEx;
stats.pEx                       = pEx;
stats.pExJit                    = pExJit;
stats.ccIn                      = ccIn;
stats.pIn                       = pIn;

stats.ttExF                     = ttExF;
stats.ttInF                     = ttInF;
stats.ccExF                     = ccExF;
stats.pExF                      = pExF;
stats.pExJitF                   = pExJitF;
stats.ccInF                     = ccInF;
stats.pInF                      = pInF;

stats.ttrelF                    = ttrelF;
stats.RprofF                    = RprofF;

details                         = stats;
details.frq_CV                  = frq;
details.cohs_CV              	= cohs;
details.phs_CV                  = phs;

stats                           = [];
stats.Q                         = [ Q Qpval ];                              % mean (over trials), SEM, p-value - geo-mean over p-values, each by spike time shuffling
stats.R                         = R;                                        % mean (over pairs), SEM, pvalue - one-sided signed rank vs. zero
stats.P                         = P;                                        % by fixed filter
stats.T                         = T;
stats.Rprof                     = Rprof;
stats.Rp                        = Rp;
stats.P2                        = P2;                                       % by new (cross-validated) filter for jitter
stats.MI                        = MI;
stats.Coh                       = mcoh;
stats.coderType                 = codeType; %0- nan, 1 - temporal, -1 - rate

return % stAnalysis

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
% h = local_calc_entropy( p )
% entropy (base 2).
%------------------------------------------------------------------------
function h = local_calc_entropy( p )

p( p == 0 )                     = NaN;
h                               = -nansum( p .* log2( p ) );

return % local_calc_entropy

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
% phi = local_circ_mean( t )
% compute mean direction
%------------------------------------------------------------------------
function phi = local_circ_mean( t )

% trigonometric functions
nans                            = isnan( t );
n                               = sum( ~nans, 1 );
x                               = cos( t );
y                               = sin( t );

% compute direction
sumx                            = nansum( x, 1 );
sumy                            = nansum( y, 1 );
C                               = sumx ./ n;
S                               = sumy ./ n;
phi                             = mod( atan2( S, C ), 2 * pi );

return % local_circ_mean

% EOF
