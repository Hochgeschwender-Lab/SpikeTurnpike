% MUTINF            mutual information from empirical distributions/counts.
%
% call              I = mutinf( mat )
%
% gets              mat     matrix of joint pdf of X, Y
%                           each element is P(x,y) or counts
%
% returns           I           p(x,y) * log2( p(x,y) / ( p(x) * p(y) ) )
%                           (sum over x,y)
%
% calls             nothing.

% 10-apr-03 ES

% last update
% August 2004

function I = mutinf( mat )

p                 	= mat / sum( sum( mat ) );                           	% p( X, Y )
px                	= sum( p, 1 );                                          % p( X )
py               	= sum( p, 2 );                                      	% p( Y )
den              	= py * px;                                              % p( X ) * p( Y )
nzi             	= p ~= 0;                                               % avoid / 0
I                   = sum( p( nzi ) .* log2( p( nzi ) ./ den( nzi ) ) );  	% actual computation

return