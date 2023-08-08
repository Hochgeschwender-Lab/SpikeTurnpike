% sta           compute sta from spike train and analog vector
%
% super simple, just mean and SD
% no error checking
%
% call                          [ m s n ev ] = sta( t, x, win )
% 
% gets                          t         event times
%                               x         analog signal
%                               win       two-element vector
%
% returns                       m         mean
%                               s         SD
%                               n         number of events in average (divide s by sqrt(n) for SEM)
%                               ev        OPTIONAL: the events used. The edge events get NaNs

% 09-jun-14 ES

% last update
% March 2015

function [ m, s, j, ev ]        = sta( t, x, win )

nt                              = length( t );                              % number of digital events
nx                              = length( x );                              % number of elements in x
nwin                            = win( 2 ) - win( 1 ) + 1;                  % number of elements in each segment

m                               = zeros( nwin, 1 );
ss                              = zeros( nwin, 1 );
if nargout > 3
    ev                          = NaN * ones( nwin, nt );
end

j                               = 0;
for i                           = 1 : nt
    t0                          = t( i ) + win( 1 );
    t1                          = t( i ) + win( 2 );
    if t0 < 1 || t1 > nx
        if nargout > 3
            if t0 < 1
                ev( -t0 + 1 + ( 1 : t1 ), i ) = x( 1 : t1 );
            elseif t1 > nx
                ev( 1 : nx - t0 + 1, i ) = x( t0 : nx );
            end
        end
        continue
    end
    j                           = j + 1;
    xi                          = x( t0 : t1 );
    if nargout > 3
        ev( :, i )              = xi;
    end
    m                           = m + xi;
    ss                          = ss + xi.^2;
end
m                               = m / j;
s                               = sqrt( ss / j - m .^ 2 );

return

% EOF