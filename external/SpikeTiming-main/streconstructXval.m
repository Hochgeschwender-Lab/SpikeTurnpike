% streconstructXval         cross-validation wrapper for streconstruct.m
%
% call                      [ yhat, cc, pval, models ] = streconstructXval( st, wn, sthat )
%
% gets                      x, x2               spike trains (see streconstruct.m)
%                           y               	analog signal
% 
% optional arguments (given as name/value pairs):
%
%                           nfolds            {10}
%                           nrepsCV           {10}
%                           minTrials         {4}; also corrected for nfolds and passed to streconstruct.m
%                           minSpikesPerTrial {1}; passed to streconstruct.m as is
%                           all other argument pairs are passed to streconstruct.m as is
% 
% returns:                  yhat              reconstruction
%                           cc                quality (mean over trials)
%                           pval              pvalue (geometric mean over trials, obtained by spike time shuffling)
%                           models            detailed information per trial
%
% calls:                    stahat, stmix, streconstruct
%                           isisort
%                           calc_spearman, mixmat, nangeomean

% 31-oct-14 ES

% last update
% August 2022

function [ yhat, cc, pval, models ] = streconstructXval( st, wn, sthat, varargin )

% initializations
yhat                            = [];
cc                              = [ NaN NaN ];
pval                            = NaN;
models                          = [];

% default values
DEFAULT_nfolds                  = 10;                                       % maximal number of folds
DEFAULT_nreps                   = 10;                                       % number of shuffles
DEFAULT_minTrials               = 4;                                        % minimal number of trials
DEFAULT_minSpikesPerTrial       = 1;
DEFAULT_verbose                 = 1;

% arguments
nargs = nargin;
if nargs < 2 || isempty( st ) || isempty( wn )
    return
end
[ m, n ]                        = size( st );
[ mw, nw ]                      = size( wn );
if m ~= mw
    error( 'input size mismatch: st and wn should have same number of rows' )
end
if n ~= nw
    if nw == 1
        wn                      = repmat( wn, [ 1 n ] );
    else
        error( 'input size mismatch: wn should be same size as st or have one column' )
    end
end
if nargs < 3
    sthat                       = [];
end
if ~isempty( sthat )
    if ~isequal( size( sthat ), [ m n ] )
        error( 'input size mismatch: st and sthat must be of same size' )
    end
end
nfolds                          = DEFAULT_nfolds;
minTrials                       = DEFAULT_minTrials;
nreps                           = DEFAULT_nreps;
minSpikesPerTrial               = DEFAULT_minSpikesPerTrial;
verbose                         = DEFAULT_verbose;
nvars                           = length( varargin );
for i                           = 1 : nvars
    if strcmpi( varargin{ i }, 'nfolds' ) && i < nvars
        nfolds                  = varargin{ i + 1 };
    end
    if strcmpi( varargin{ i }, 'mintrials' ) && i < nvars
        minTrials               = varargin{ i + 1 };
    end
    if strcmpi( varargin{ i }, 'nrepsCV' ) && i < nvars
        nreps                   = varargin{ i + 1 };
    end
    if strcmpi( varargin{ i }, 'minSpikesPerTrial' ) && i < nvars
        minSpikesPerTrial       = varargin{ i + 1 };
    end
    if strcmpi( varargin{ i }, 'verbose' ) && i < nvars
        verbose                 = varargin{ i + 1 };
    end
end

if nfolds == 0
    return
end

% partition into train and test trials
nfolds                          = min( nfolds, n );
minTrialsHat                    = minTrials - ceil( minTrials / nfolds );
midx                            = mixmat( ( 1 : n )', 1, 1 )';

% prepare for streconstruct.m call
ridx                            = [];
for i                           = 1 : nvars
    if strcmpi( varargin{ i }, 'nfolds' ) && i < nvars
        ridx                    = [ ridx i + [ 0 1 ] ];
    end
    if strcmpi( varargin{ i }, 'nrepsCV' ) && i < nvars
        ridx                    = [ ridx i + [ 0 1 ] ];
    end
    if strcmpi( varargin{ i }, 'mintrials' ) && i < nvars
        varargin{ i + 1 }       = minTrialsHat;
    end
    if strcmpi( varargin{ i }, 'verbose' ) && i < nvars
        varargin{ i + 1 }       = 0;
    end
end
varargin( ridx )                = [];

yhat                            = NaN * ones( m, n );
cc                              = NaN * ones( 1, n );
pvals                           = NaN * ones( 1, n );
model1                          = struct( 'f', [], 'isiedges', [], 'optlag', [] ...
    , 'trainset', [], 'testset', [], 'pval', [], 'cc', [] );
models                          = repmat( model1, [ nfolds 1 ] );

if verbose
    fprintf( '%s: %d-fold cross-validation (%d randomizations/trial) ', upper( mfilename ), nfolds, nreps )
end

% cross validate:
% -build model for training set
% -apply for test set
% -shuffle the test set nrepsCV times, and re-apply the cross-validated model

for k                           = 1 : nfolds
    
    % indices (training/test trials)
    testidx                     = sort( midx( k : nfolds : n ) );           % randomized
    trainidx                    = setdiff( 1 : n, testidx );                % training set
    if all( full( sum( st( :, testidx ) ) ) < minSpikesPerTrial )
        continue
    end
    if verbose
        fprintf( '.' )
    end
    
    % training set (nt-1, or floor(nt/2) trials)
    if isempty( sthat )
        st2                     = sthat;
    else
        st2                     = sthat( :, trainidx );
    end
    [ ~, statsK ]               = streconstruct( st( :, trainidx ), wn( :, trainidx ), 'x2', st2, varargin{ : } );
    
    % model
    if isempty( statsK )
        continue
    end
    f                           = statsK.f;
    isiedges                    = statsK.isiedges;
    optlag                      = statsK.optlag;
    
    % test set (xk) + shuffles (xm)
    xk                          = st( :, testidx );
    xm                          = stmix( xk, [], nreps );
    xt                          = isisort( [ xk xm ], isiedges );
    yhatK                       = stahat( xt, f, m, optlag );
    yhat( :, testidx )          = yhatK( :, 1 : length( testidx ) );
    
    % reconstruction quality (per trial)
    wnk                         = wn( :, testidx );
    wnidx                       = ones( nreps, 1 ) * ( 1 : length( testidx ) ) ;
    wnm                       	= wnk( :, wnidx( : ) );
    wnK                         = [ wnk wnm ];
    
    rCV                         = calc_spearman( yhatK, wnK );
    rRand                       = rCV( ( length( testidx ) + 1 ) : end );
    cc( testidx )               = rCV( 1 : length( testidx ) );
    
    % reconstruction p-value (per trial)
    for ki                      = 1 : length( testidx )
        cidx                    = ( 1 : nreps ) + ( ki - 1 ) * nreps;
        if nreps < 100
            mm                  = nanmean( rRand( :, cidx ) );
            ss                  = nanstd( rRand( :, cidx ) );
            pvals( testidx( ki ) )    = 1 - normcdf( cc( testidx( ki ) ), mm, ss );
        else
            pvals( testidx( ki ) )    = ( sum( rRand >= cc( testidx( ki ) ) ) + 1 ) / ( nreps + 1 );
        end
    end
    
    % summarize
    models( k )                 = struct( 'f', f, 'isiedges', isiedges, 'optlag', optlag ...
        , 'trainset', trainidx, 'testset', testidx, 'pval', pvals( testidx ), 'cc', cc( testidx ) );

end

pval                            = nangeomean( pvals + eps );
cc                              = [ nanmean( cc, 2 ) std( cc, [], 2 ) / sqrt( size( cc, 2 ) ) ];
if verbose
    fprintf( 'cc=%0.3g; p=%0.3g\n', cc( 1 ), pval )
end

return

% EOF
