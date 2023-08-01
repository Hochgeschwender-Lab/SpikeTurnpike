% streliability           	estimate reliability of single-cell spike trains
%
% call                      [ rel, ccSD, ccSDF, sF ] = streliability( x, yhat )
%
% gets                      x               multi-trial spike trains (sparse matrix, n x nt)
%                           yhat            a specific reconstruction (same dimensions as x)
% 
% optional arguments (given as name/value pairs):
%
%                           s               {[ 1 2 4 8 16 32 64 ]}  m-element vector of Gaussian SDs [samples]
%                           ftype           { 'gauss' }             filter type; 'exp', 'alpha', and 'rect' are also supported
%                           doFocused       {1}                     flag
%                           graphics        {0}                     flag
%
% returns                   rel             mean, SEM, and p-value of all cc (for all possible yhat pairs); p-value by one-sided signed-rank test
%                           ccSD            for all pairs, and all jitters (m by n*(n-1)/2 matrix)
%                           ccSDF           for all pairs, and all jitters, for a focused range of SDs
%                           sF              the Gaussian SDs used for the focused range
%
% calls                     ParseArgPairs
%                           calc_spearman
%
% does                      answers the following two questions:
%                               1. how reliable is the reconstruction yhat
%                               2. how reliable are the spike trains x (in terms of reconstruction with a non data-dependent filter)
%
% definition:               reliability := reproducibility of the same spike train
%                                   this is similar to the Van Rossum metric (Neural Computation 2001) 
%                                   and identical to the Schreiber et al. metric (Neurocomputing, 2003).
%

% 27-oct-14 ES

% last update
% August 2022

function [ rel, ccSD, ccSDF, sF ] = streliability( x, yhat, varargin )

% block-wise to prevent memory exhaustion:
BLOCKSIZE                       = 50;

% arguments
nargs                           = nargin;
if nargs < 2 || isempty( x ) || isempty( yhat )
    return
end
[ s, ftype ...
    , doFocused ...
    , graphics ]                = ParseArgPairs(...
    { 's', 'ftype' ...
    , 'doFocused' ...
    , 'graphics' }...
    , { 2 .^ ( 0 : 6 ), 'gauss' ...
    , 1 ...
    , 0 }...
    , varargin{ : } );
nt                              = size( x, 2 );
if ~isequal( size( x ), size( yhat ) )
    error( 'input size mismatch' )
end

% preps
aidx                            = local_make_ai( nt );
npairs                          = size( aidx, 1 );
blocks                          = local_makeblocks( size( aidx, 1 ), BLOCKSIZE );
nblocks                         = size( blocks, 1 );

%------------------------------------------------------------------------
% (1) compute reliability for the reconstruction:
cc                              = zeros( npairs, 1 );
for i                           = 1 : nblocks
    bidx                        = blocks( i, 1 ) : blocks( i, 2 );
    cc( bidx )                  = calc_spearman( yhat( :, aidx( bidx, 1 ) ), yhat( :, aidx( bidx, 2  ) ) );
end

% compute p-value
relPval                         = signrank( cc, 0 );
if nanmedian( cc ) > 0
    relPval                     = relPval / 2;
else
    relPval                     = 1 - relPval / 2;
end
rel                             = [ nanmean( cc ) std( cc ) / sqrt( size( cc, 1 ) ) relPval ];

%------------------------------------------------------------------------
% (2) compute reliability for various temporal jitters:
nSD                             = length( s );
if all( isnan( s ) )
    nSD                         = 0;
end
ccSD                            = NaN * ones( nSD, npairs );
xf                              = full( x );

for k                           = 1 : nSD

    % build the filter
    tau                         = s( k );
    switch lower( ftype )
        case 'gauss'
            iwin                = local_makegaussfir( tau );
        case 'alpha'
            t                   = 0 : 1 : 5 * tau;
            iwin                = t .* exp( -t / tau );
        case 'exp'
            t                   = 0 : 1 : 5 * tau; % 3 is >95% of the support, 5 is >99% of the support
            iwin                = exp( -t / tau );
        case { 'rect', 'ma', 'box', 'boxcar' }
            iwin                = ones( tau, 1 );
        otherwise
            error( 'unsupported filter type' )
    end
    iwin                        = iwin / sum( iwin );
    
    % convolve the spike trains:
    switch lower( ftype )
        case { 'exp', 'alpha' }
            fg                  = filter( iwin, 1, xf );                    % exponential tail
        case { 'gauss', 'rect', 'ma', 'box', 'boxcar' }
            fg                  = local_firfilt( xf, iwin );                % symmetric around spike
    end
    
    % compute the cc:
    for i                       = 1 : nblocks
        bidx                    = blocks( i, 1 ) : blocks( i, 2 );
        ccSD( k, bidx )         = calc_spearman( fg( :, aidx( bidx, 1 ) ), fg( :, aidx( bidx, 2  ) ) );
    end
    
end

%------------------------------------------------------------------------
% (3) iterate around peak of R profile:
Rprof                           = mean( ccSD, 2 );
[ ~, snum ]                     = max( Rprof );

if s( snum ) > 1 && doFocused
    
    % repeat everything for the range between s( maxidx ) and the one before it
    sF                          = ( s( snum - 1 ) : s( snum ) )';
    nFocused                    = length( sF );
    ccSDF                       = NaN * ones( nFocused, npairs );
    
    for k                       = 1 : length( sF )
        
        % build the filter
        tau                     = sF( k );
        switch lower( ftype )
            case 'gauss'
                iwin            = local_makegaussfir( tau );
            case 'alpha'
                t               = 0 : 1 : 5 * tau;
                iwin            = t .* exp( -t / tau );
            case 'exp'
                t               = 0 : 1 : 5 * tau;                      	% 3 is >95% of the support, 5 is >99% of the support
                iwin            = exp( -t / tau );
            case { 'rect', 'ma', 'box', 'boxcar' }
                iwin            = ones( tau, 1 );
            otherwise
                error( 'unsupported filter type' )
        end
        iwin                    = iwin / sum( iwin );
        
        % convolve the spike trains:
        switch lower( ftype )
            case { 'exp', 'alpha' }
                fg              = filter( iwin, 1, xf );                    % exponential tail
            case { 'gauss', 'rect', 'ma', 'box', 'boxcar' }
                fg              = local_firfilt( xf, iwin );              	% symmetric around spike
        end
        
        % compute the cc:
        for i                   = 1 : nblocks
            bidx                = blocks( i, 1 ) : blocks( i, 2 );
            ccSDF( k, bidx )    = calc_spearman( fg( :, aidx( bidx, 1 ) ), fg( :, aidx( bidx, 2  ) ) );
        end
        
    end
    
else
    
    sF                          = [];
    ccSDF                       = [];
    
end


%------------------------------------------------------------------------
% (4) plot if required
if ~graphics
    return
end

newplot
Rsem                            = std( ccSD, [], 2 ) / sqrt( size( ccSD, 2 ) );
if nSD ~= 0
    plot( log2( s ), Rprof, 'b', log2( s ), Rprof + Rsem, '--b' ...
        , log2( s ), Rprof - Rsem, '--b' );
    xlim( [ min( log2( s ) ) max( log2( s ) ) ] )
    xlims                       = [ min( log2( s ) ) max( log2( s ) ) ];
    ticks                       = floor( xlims( 1 ) ) : ceil( xlims( 2 ) );
    set( gca, 'xlim', xlims, 'xtick', ticks, 'xticklabel', 2.^ticks  )
    xlabel( 's [samples]' )
end
ylabel( 'Reliability (cc)' )
line( xlim, rel( 1 ) * [ 1 1 ], 'color', [ 1 0 0 ], 'linestyle', '-' )
line( xlim, rel( 1 ) + [ -1 -1 ] * rel( 2 ), 'color', [ 1 0 0 ], 'linestyle', '--' )
line( xlim, rel( 1 ) + [ 1 1 ] * rel( 2 ), 'color', [ 1 0 0 ], 'linestyle', '--' )
axis tight
set( gca, 'tickdir', 'out', 'box', 'off' ),

return % streliability

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
% K = local_makegaussfir( sigmaX, N )
% create a 1D gaussian FIR
%------------------------------------------------------------------------
function K = local_makegaussfir( sigmaX )

n                               = 2 * floor( sigmaX / 2 ) + 1; % make odd
support                         = 6;
N                               = support * n + 1;
x                               = -( N - 1 ) / 2 : ( N - 1 ) / 2;
K                               = 1 / ( 2 * pi * sigmaX ) * exp( -( x.^2 / 2 / sigmaX^2 ) );
K                               = K / sum( K( : ) );

return % local_makegaussfir

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

% EOF
