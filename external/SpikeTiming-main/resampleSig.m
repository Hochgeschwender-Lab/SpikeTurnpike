% resampleSig       wrapper for resample.m
%
% does:             works with sampling rates and takes care of edges
%
% call:             y = resampleSig( x, Fs, newFs )
%
% gets:             x: signal (treated as a vector)
%                   Fs: sampling rate [samples/s]
%                   newFs: desired sampling rate [samples/s]
%
% returns:          y: ceil( length( x ) * newFs / Fs ) elements
%
% use example:
% >> Fs = 10000; T = 1; WN = randn( T * Fs, 1 );
% % resample to 25000:
% >> y = resampleSig( WN, Fs, 25000 )
%
% to plot:
% >> figure, plot( 1 / Fs : 1 / Fs : T, WN, '.-b', 1 / Fs : 1 / newFs : ( T + 1 / Fs - 1 / newFs ) , y, '-r' )

% 22-mar-14 ES

% last update
% August 2022

function y = resampleSig( x, Fs, newFs )

y                               = [];

nargs = nargin;
if nargs < 1 || isempty( x )
    return
end
if nargs < 2 || isempty( Fs )
    Fs                          = 1;
end
if nargs < 3 || isempty( newFs )
    newFs                       = 1;
end

% get the rational numbers
str                             = rats( newFs/Fs );
nstr                            = length( str );
tok                             = strfind( str, '/' );
if isempty( tok )
    pretok                      = str;
    Q                           = 1;
else
    pretok                      = str( 1 : ( tok - 1 ) );
    posttok                     = str( ( tok + 1 ) : nstr );
    Q                           = str2double( posttok );
end
P                               = str2double( pretok );

% work with a column vector
x                               = x( : );

% mirror-expand
nx                              = length( x );
xr                              = x( nx : -1 : 1 );
xpad                            = [ xr; x; xr ];

% resample
ypad                            = resample( xpad, P, Q );

% clip back
if ( length( ypad ) / 3 ) ~= ceil( nx * P / Q )
    %warning( 'resample assumptions violated - check!!' )
end
ny                              = round( length( ypad ) / 3 );
y                               = ypad( ( ny + 1 ) : ( 2 * ny ) );

return

% EOF