% stahat                        reconstruct analog signal
%
% call                          xhat = stahat( t, f, n, lag )
% 
% gets                          t           sparse matrix (multiple columns). individual elements can be tagged
%                               f           filter, can be a filter bank
%                               n           number of elements in output
%                               lag         to shift the output
%
% see also:                     isisort, sta, stacalc

% 27-oct-14 ES

% last update
% August 2022

function xhat                   = stahat( t, f, n, lag )

% arguments
nargs                           = nargin;
if nargs < 2 || ~issparse( t ) || isempty( f )
    error( 'input mismatch' )
end
nt                              = size( t, 2 );
nbins                           = size( f, 2 );
if nargs < 3 || isempty( n )
    n                           = max( t ) + ceil( size( f, 1 ) / 2 );
end
if nargs < 4 || isempty( lag )
    lag                         = 0;
end

% initialize
if nbins > 1
    utags                       = setdiff( full( unique( t ) ), 0 );
    if ~isempty( setdiff( utags, 0 : nbins ) )
        error( 'filter/data mismatch' )
    end
    ntags                       = length( utags );
end
xhat                            = zeros( n, nt );

% reconstruct
if nbins == 1
    xhat                        = local_firfilt( full( t ), f );            % enables t(i,j)>1
else
    for i                       = utags( : ).'
        idx                     = t == i;
        xi                      = sparse( n, nt );
        xi( idx )               = 1;
        xhat                    = xhat + local_firfilt( full( xi ), f( :, i ) );
    end
end

% shift 
if lag ~= 0
    xhat                        = local_shift( xhat, lag, 0 );
end

return % stahat

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

% EOF
