% calc_fwhh             determine full-width at half-height (direct method)
%
% call                  [ width, lag, val ] = calc_fwhh( x )
%
% gets                  x       signal. if x is a matrix, works on columns
% 
% returns               width       FWHH [samples]
%                       lag         of extremum [samples]
%                       val         of the extremum
% 
% NOTE:                 assumes a zero reference
% 
% calls                 nothing

% 23-mar-14 ES

% last update
% April 2020

function [ width, lag, val ] = calc_fwhh( x )

% arguments
if nargin < 1 || isempty( x )
    return
end
[ m, n ]                            = size( x );
width                               = NaN * ones( 1, n );
lag                                 = width;
val                                 = width;
idx                                 = NaN * ones( m, n );
if all( isnan( x ) )
    return
end

% computations
ref                                 = zeros( 1, n );
midx                                = max( abs( x ), [], 1 ) == max( x, [], 1 );
if sum( midx ) > 0
    [ val( midx ), lag( midx ) ]    = max( x( :, midx ) );
    halfAmp                         = ( val( midx ) + ref( midx ) ) / 2;
    idx( :, midx )                  = cumsum( bsxfun( @lt, x( :, midx ), halfAmp ) );
end
midx                                = max( abs( x ), [], 1 ) == -min( x, [], 1 );
if sum( midx ) > 0
    [ val( midx ), lag( midx ) ]    = min( x( :, midx ) );
    halfAmp                         = ( val( midx ) + ref( midx ) ) / 2;
    idx( :, midx )                  = cumsum( bsxfun( @gt, x( :, midx ), halfAmp ) );
end

nnans                               = ~isnan( lag );
idxvalid                            = idx( :, nnans );
eidx                                = idxvalid( lag( nnans ) + m * ( 0 : ( sum( nnans ) - 1 ) ) );
width( nnans )                      = sum( bsxfun( @eq, idxvalid, eidx ) ) - 1;

return

% EOF

