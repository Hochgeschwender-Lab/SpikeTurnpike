% stwiener              construct a modified Wiener filter for spike train and analog signal
%
% call                  [ f fsem ] = stwiener( t, x, win )
% 
% gets                  t                   spike times
%                       x                   analog signal
%                       win                 window
%
% returns               f                   filter: R(x,t)/R(t,t)
%                       fsem                filter SEM 
%
% does                  computes the Z-scored modified Wiener filter and its SEM
%
% note:
% should pre-process the data so that there are no spikes at a distance
%       of max( abs( win ) ) from the edges of x
%
% calls                 sta
%                       calc_cor (mex)

% 14-jun-15 ES

% last update
% August 2022

function [ f, fsem, f0 ] = stwiener( t, x, win )

% arguments
nargs                           = nargin;
if nargs < 3 || isempty( t ) || isempty( x ) || isempty( win )
    error( 'input missing' )
end
if issparse( t ) || ~isequal( t, t( : ) )
    error( 't must be a list of spike times' )
end
if isempty( x ) || ~isequal( x, x( : ) )
    error( 'x must be an analog vector' )
end
switch length( win )
    case 1
        win                     = [ -abs( win ) abs( win ) ];
    case 2
        if ~isequal( win, sort( win ) )
            error( 'win must be sorted' )
        end
    otherwise
        error( 'win must be a two-vector element' )
end

% preps
maxlag                          = max( abs( win ) );

% computation
[ Rtx, ss, nn ]                 = sta( t, x( : ), win );                 	% input-output cross-corr
rtt                             = calc_cor( t, t, maxlag );               	% input auto-corr
tidx                            = [ maxlag + ( 1 : maxlag + 1 ) 1 : maxlag ];       % rearrange (non-causal, two sided)
[ ~, i1 ]                       = intersect( win( 1 ) : win( 2 ), -maxlag : maxlag );
tidx                            = tidx( i1 );                             	% support causal windows
Rtt                             = toeplitz( rtt( tidx ) );                  % input correlation matrix
f0                              = Rtt \ Rtx;                              	% solve the equation

% filter SEM
fsem0                           = ss / sqrt( nn ) / max( rtt );

% scale both to Z-scores
mf                              = mean( f0 );
sf                              = std( f0 );
f                               = ( f0 - mf ) / sf;
fsem                            = ( fsem0 - mf ) / sf;

return

% EOF