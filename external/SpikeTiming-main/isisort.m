% isisort               sort spikes into ISI groups
%
% call                  [ xc bedges cpb ] = isisort( x, bedges, x2 )
% 
% gets                  x           spike train/s. can be either a sparse matrix or a vector of spike times
%                       n           number of ISI bins (scalar) or an n+1 vector of ISI edges
%                       x2          (optional) reference spike train (cross-cell model, e.g. xsta or xwiener)
%
% returns               xc          classified spike trains, either a sparse matrix of a vector
%                                       (depending on the input format)
%                       bedges      edges; (n+1 x 1 vector)
%                       cpb         summary of counts per bin (n x 1 vector)
%
% calls                 nothing
%
% see also              streconstruct

% 27-oct-14 ES

% last update
% August 2022

function [ xc, bedges, cpb ] = isisort( x, bedges, x2 )

nargs                           = nargin;
if nargs < 1 || isempty( x ) || issparse( x ) && sum( x( : ) ) == 0
    xc                          = x;
    bedges                      = [];
    cpb                         = [];
    return
end
if nargs < 2 || isempty( bedges )
    bedges                      = 1;
end
if nargs < 3
    x2                          = [];
end

ftype                           = '';
if isempty( x2 )
    if issparse( x ) || isvector( x )
        ftype                   = 'cmodel';
    end
else
    if issparse( x ) && issparse( x2 ) && isequal( size( x ), size( x2 ) ) ...
            || ~issparse( x ) && ~issparse( x2 ) && isvector( x ) && isvector( x2 )
        ftype                   = 'xmodel';
    end
end

% 1. calc the ISIs
if isequal( ftype, 'cmodel' )
    isis                        = local_calc_isis( x );
else
    isis                        = local_calc_isis( x, x2 );
end

% 2. get the ISI edges
if isempty( bedges ) || length( bedges ) == 1 && bedges( 1 ) < 2
    xc                          = x;
    cpb                         = [];
    return
elseif length( bedges ) == 1
    nbins                       = bedges;
    [ sval, eval ]              = local_make_equal_bins( isis, nbins );
    nbins                       = length( sval );
    eval( nbins )               = inf;
    bedges                      = [ sval; eval( nbins ) ];
else
    if any( diff( round( bedges ) ) <= 0 )
        error( 'bedges must be a vector of increasing integers' )
    end
    nbins                       = length( bedges ) - 1;
end
isis( isnan( isis ) )           = inf;

% 3. sort the spikes (first spike in trial gets inf)
[ cpb, sidx ]                   = histc( isis, bedges );
cpb( nbins )                    = sum( cpb( end + [ -1 0 ] ) );
cpb( nbins + 1 )                = [];
sidx( sidx == nbins + 1 )       = nbins;

if issparse( x )
    xc                          = x;
    xc( find( x( : ) ) )        = sidx;                                     % keep the classified trains
else
    xc                          = sidx;
end

return % isisort

%------------------------------------------------------------------------
% y = local_calc_isis( st1, st2, pad )
% or cross-isis
%------------------------------------------------------------------------
function y = local_calc_isis( st1, st2, pad )

y                               = [];
nargs                           = nargin;
if nargs < 2
    st2                         = [];
end
if nargs < 3 || isempty( pad )
    pad                         = NaN;
end
pad                             = pad( 1 );

if issparse( st1 )
    ntrials                     = size( st1, 2 );
    for i                       = 1 : ntrials
        x1                      = find( st1( :, i ) );
        if isempty( x1 )
            continue
        end
        if ~isempty( st2 ) && issparse( st1 ) && isequal( size( st1 ), size( st2 ) )
            x2                  = find( st2( :, i ) );
        else
            x2                  = [];
        end
        y1                      = calc_isis1( x1, x2, pad );
        y                       = [ y; y1 ];
    end
else
    st1                         = st1( : );
    if isempty( st1 ) 
        return
    end
    if ~isempty( st2 )
        st2                     = st2( : );
    end
    y                           = calc_isis1( st1, st2, pad );
end
                
return

%------------------------------------------------------------------------
% y = calc_isis1( st1, st2, pad )
% callback from calc_isis
%------------------------------------------------------------------------
function y = calc_isis1( st1, st2, pad )

nargs                           = nargin;
if nargs < 2
    st2                         = [];
end
if ~isempty( st1 ) && ~issorted( st1 )
    st1                         = sort( st1 );
end
if isempty( st2 )
    y                           = [ pad; diff( st1 ) ];
    return
end
if ~issorted( st2 )
    st2                         = sort( st2 );
end
i2                              = 1;
y                               = pad * ones( size( st1 ) );
for i1                          = 1 : length( st1 )
    s1                          = st1( i1 );
    if i2 < length( st2 ) && s1 > st2( i2 + 1 )
        i2                      = i2 + 1;
    end
    s2                          = st2( i2 );
    if s1 >= s2                                                             % synchrony allowed
        y( i1 )                 = s1 - s2;
    end
end

return % calc_isis1

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

% EOF
