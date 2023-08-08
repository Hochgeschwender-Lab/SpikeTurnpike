% stmix                 shuffle/jitter spikes in a sparse array
%
% call                  xm = stmix( x, halfwin, nreps, jmode )
%
% gets                  x               (m x n) sparse array
%                       halfwin         {[]}; [samples]. empty shuffles
%                       nreps           {1}; if an integer > 1, replicates each column nreps times
%                       jmode           {1}     0:      shuffle (halfwin ignored)
%                                               1:      spike time jitter (-h to h)
%                                               2:      interval jitter (2*h+1 interval)
%
% returns               xm              m x ( n * nreps ) randomized sparse array, 
%                                           organized in blocks of nreps columns:
%                                       first nreps columns of xm are randomized versions of the 
%                                       first column of x
%
% calls                 nothing
%
% see also              streconstruct

% 17-jun-14 ES

% last update
% October 2020

function xm = stmix( x, h, p, j, s )

xm                              = [];

% arguments
nargs                           = nargin;
if nargs < 1 || isempty( x )
    return
end
if nargs < 2 || isempty( h )
    h                           = [];
end
if nargs < 3 || isempty( p )
    p                           = 1;
end
if nargs < 4 || isempty( j )
    j                           = 1;
end
if nargs < 5 || isempty( s )
    s                           = 1;
end

if isempty( h ) 
    j                           = 0;
end
if ~isempty( h ) && ( h < 0 || isinf( h ) )
    error( 'halfwin must be non-negative number' )
end
if p < 0 || p ~= round( p ) || isinf( p )
    error( 'nreps must be a non-negative integer' )
end
if ~ismember( j, 0 : 2 )
    error( 'j must be 0, 1, or 2' )
end
if ~ismember( s, 0 : 1 )
    error( 'j must be 0 or 1' )
end

% intialize
if ~issparse( x )
    return
end
[ m, n ]                        = size( x );
nspks                           = full( sum( x ) );
cspks                           = [ 0 cumsum( p * nspks ) ];
ntot                            = sum( nspks );
trN                             = zeros( ntot * p, 1 );
if s
    rng( round( rand( 1 ) * sum( 100 * clock ) ), 'v4' )
end

% randomize
switch j
    case 0                                                                  % shuffle (re-distribute over duration)
        [ ~, rnd ]              = sort( rand( m, n * p ), 1 );
        for i                   = 1 : n
            cols                = ( 1 : p ) + ( i - 1 ) * p;
            idx                 = ( cspks( i ) + 1  ) : cspks( i + 1 );
            d1                  = ones( nspks( i ), 1 ) * ( cols - 1 ) * m;
            tr                  = rnd( 1 : nspks( i ), cols ) + d1;
            trN( idx )          = tr( : );
        end
        
    case 1                                                                  % jitter (spike-time)
        jit                     = round( 2 * h * rand( ntot * p, 1 ) - h );
        for i                   = 1 : n
            cols                = ( 1 : p ) + ( i - 1 ) * p;
            idx                 = ( cspks( i ) + 1  ) : cspks( i + 1 );
            d1                  = ones( nspks( i ), 1 ) * ( cols - 1 ) * m;
            st                  = find( x( :, i ) ) * ones( 1, p );
            d                   = jit( idx );
            trN( idx )          = st( : ) + d1( : ) + d;
        end
        
    case 2                                                                  % jitter (interval)
        [ row, col ]            = find( x );
        st                      = sub2ind( [ m n ], row, col );
        wn                      = ceil( st / ( 2 * h + 1 ) );
        r                       = ceil( rand( ntot, p ) * ( 2 * h + 1 ) );
        nst                     = ( repmat( wn, [ 1 p ] ) - 1 ) * ( 2 * h + 1 ) + r;
        m1                      = ones( ntot, 1 ) * ( 0 : 1 : ( p - 1 ) ) * m;
        m2                      = ( col - 1 ) * m * ( p - 1 ) * ones( 1, p );
        trN                     = m1 + m2 + nst;
        
end

% build arrays
trN( trN < 1 )                  = [];
trN( trN > ( m * n * p ) )      = [];
[ ii, jj ]                      = ind2sub( [ m n * p ], trN );
xm                              = sparse( ii, jj, 1, m, n * p );

return

% EOF