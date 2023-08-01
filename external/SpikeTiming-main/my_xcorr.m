% MY_XCORR          compute the column-wise cross-correlation between real matrices.
%
% call              [ cc, lags ] = my_xcorr( x, y, maxlag, normf )
%
% gets              x           matrix, trials in columns
%                   y           {x}; another matrix of same dimensions; if
%                                   left empty, computes auto-correlation
%                   maxlag      {size(x,2)}; maximal lag
%                   normf       {0}     raw cross correlation
%                               1       scale by autocorrelation at tau 0 (Pearson correlations) 
%                               -1      cross-covariances (remove means)
%
% does              column-wise cross-correlation
%
% NOTE:             if y lags after x, peak is in negative lags
%                   (i.e. y is the REFERENCE signal)
%
% calls             nothing

% 23-feb-05 ES

% last update
% April 2020

function [ cc, lags ] = my_xcorr( x, y, maxlag, normf )

% arguments
nargs                           = nargin;
if nargs < 1 || isempty( x )
    error( 'need input!' )
end
if nargs < 2 || isempty( y )
    y                           = x;
end
if nargs < 3
    maxlag                      = []; 
end
if nargs < 4
    normf                       = 0; 
end
if ~isequal( size( x ), size( y ) )
    error( 'X, Y must be equal size' );
end
if ~isempty( maxlag )
    if maxlag < 0 || maxlag ~= round( maxlag )
        error( 'MLAG must be non-negative scalar' ); 
    end
end

% preparation
n                               = size( x, 1 );
if isempty( maxlag )
    maxlag                      = n; 
end
MLAG0                           = maxlag;
lags                            = -maxlag : maxlag;
maxlag                          = min( n - 1, maxlag );
if normf < 0                                                                % compute cross-covariances
    x                           = x - ones( n, 1 ) * mean( x );
    y                           = y - ones( n, 1 ) * mean( y );
end

% computation
nfft                            = 2 ^ nextpow2( 2 * n - 1 );
xx                              = fft( x, nfft );
yy                              = fft( y, nfft );
cc                              = real( ifft( xx .* conj( yy ) ) );
cc                              = [ cc( nfft - maxlag + 1 : nfft, : ); cc( 1 : maxlag + 1, : ) ];

if normf % compute Pearson correlations (or cross-covariances)
    if ~isequal( x, y ) 
        % distinct signals: normalize by autocorrelations at zero lag
        cxx0                    = sum( abs( x ) .^ 2 );
        cyy0                    = sum( abs( y ) .^2 );
        scale                   = sqrt( cxx0 .* cyy0 );
        cc                      = cc ./ repmat( scale, size( cc, 1 ), 1 );
    else
        % same signal (autocorrelation case), normalize by c[0]
        cc                      = cc ./ repmat( cc( maxlag + 1, : ), size( cc, 1 ), 1 );
    end
end

if maxlag < MLAG0
    pad                         = zeros( MLAG0 - maxlag, size( cc, 2 ) );
    cc                          = [ pad; cc; pad ];
end

return

% EOF