% stcodertypes       	classify reliability profiles into rate or temporal coders
%
% call                  [ codeType, max_out ] = stcodertypes( mm, ss, nn, bb, graphics )
%
% gets                  mm          mean realibility (m x n matrix)
%                       ss          SEM of the realibility (m x n matrix)
%                       nn          number of trials (not pairs!), n-element vector
%
% optional arguments (given as a list):
%
%                       bb          {1:length(nn)}  bin values (Gaussian SDs), m-element vector
%                       rFull       {[]}            all pairs
%                       graphics    {1}             1   plots 3 panels in a new figure
%                                                   -1  plots in current axes (newplot)
%                                                   0   does not plot
%
% returns               codeType    n-element vector:   1   temporal coder
%                                                      -1   rate coder
%                                                       0   unclassified
%                       max_out     sigma in which local max of each Rprof
%
% temporal coder    := peak in curve before the last bin
%       (1) local peak (max) at any bin before the last bin
%       (2) must be sig. above adjacent two bins
%
% rate coder        := ~increasing R with delta (dR/dDelta >= 0 for all Delta)
%       (1) global max is at last bin
%       (2) max is sig. above zero

% 15-oct-19 ES

% last update
% August 2022

function [ codeType, max_out ] = stcodertypes( mm, ss, nn, bb, rFull )

% constants
pTH                             = 0.05;

% arguments
nargs                           = nargin;
if nargs < 3 || isempty( mm ) || isempty( ss ) || isempty( nn )
    return
end
[ nBins, nUnits ]               = size( mm );
if nargs < 4 || isempty( bb )
    bb                          = 1 : nBins;
end
if nargs < 5 || isempty( rFull )
    rFull                       = [];
end
if ~isequal( size( mm ), size( ss ) ) || nUnits ~= length( nn ) || nBins ~= length( bb )
    error( 'input size mismatch' )
end

% initialize output
codeType1                       = zeros( nUnits, 1 );
codeType2                       = zeros( nUnits, 1 );
max_out                         = zeros( nUnits, 1 );
max_loc                         = zeros( nUnits, 1 );

% go over profiles
for i                           = 1 : nUnits
    
    %------ global maximum
    [ ~, maxidx ]               = max( mm( :, i ) );
    max_out( i )                = bb( maxidx );
    % do a 1-sided ttest to quickly check if significantly larger than zero
    np                          = nn( i ) * ( nn( i ) - 1 ) / 2;
    [ ~, pp ]                   = local_ttesthat( mm( :, i ), ss( :, i ), np, 0, [], 1 );
    
    % if at second-to-last bin, check if sig. larger than last bin
    if pp( maxidx ) <= pTH && mm( maxidx, i ) > 0
        if maxidx < nBins                                                   % temporal coder
            OK                  = 1;
            if maxidx == ( nBins - 1 )
                if ~isempty( rFull )
                    pp2         = ranksum( rFull( maxidx, : ), rFull( maxidx + 1, : ) );
                else
                    pp2         = local_ttest2hat( mm( maxidx, i ), ss( maxidx, i ) ...
                        , nn( i ), mm( maxidx + 1, i ), ss( maxidx + 1, i ), np );
                end
                if pp2 > pTH
                    OK          = 0;
                end
            end
            if OK
                codeType1( i )  = 1;
            end
        else                                                                % rate coder
            codeType1( i )      = -1;
        end
    end
    
    %------ local maximum
    % look for local maxima distinct from the global max and not in last bin
    if ~isempty( rFull )
        mm1                     = mm;
        if ~isvector( mm )
            [ mm, Imm ]         = unique( mm', 'rows', 'stable' );
            mm                  = mm';
        else
            [ mm, Imm ]         = unique( mm, 'rows', 'stable' );
        end
    end
    cidx                        = local_find_local_extrema( mm( :, i ), 'max' );
    cidx                        = setdiff( cidx, [ maxidx nBins ] );
    [ ~, maxidxL ]              = max( mm( cidx, i ) );
    maxidxL                     = cidx( maxidxL );
    if ~isempty( maxidxL )
        % find the nearest local minimum to the right
        mins                    = local_find_local_extrema( mm( :, i ), 'min' );
        idx                     = find( mins > maxidxL, 1, 'first' );
        if ~isempty( idx )
            if ~isempty( rFull )
                minidxL1        = mins( idx );
                minidxL         = Imm( minidxL1 );
                mm              = mm1;
                pp3             = ranksum( rFull( maxidxL, : ), rFull( minidxL, : ) );
            else
                minidxL         = mins( idx );
                pp3             = local_ttest2hat( mm( maxidxL, i ), ss( maxidxL, i ) ...
                    , nn( i ), mm( minidxL, i ), ss( minidxL, i ), np );
            end
            if pp3 <= pTH
                codeType2( i )  = 1;
                max_loc( i )    = bb( maxidxL);
            end
        end
    end
    if ~isempty (rFull)
        mm                      = mm1;
    end
    
end

% add logic
codeType                        = codeType1;
codeType( codeType1 == -1 & codeType2 == 1 ) = 1;

for i                           = 1 : length( max_out )
    if max_loc( i ) > 0
        max_out( i )            = max_loc( i );
    end
end

return % stcodertypes

%------------------------------------------------------------------------
% [ h, p, ci, stats ] = local_ttesthat( mx, sx, nx, m, alpha, tail )
% identical to ttest, works on summary stats (rather than raw data)
%------------------------------------------------------------------------
function [ h, p, ci, stats ] = local_ttesthat( mx, sx, nx, m, alpha, tail )

if nargin < 4 || isempty(m)
    m                           = 0;
end

if nargin < 5 || isempty(alpha)
    alpha                       = 0.05;
elseif ~isscalar(alpha) || alpha <= 0 || alpha >= 1
    error(message('stats:ttest:BadAlpha'));
end

if nargin < 6 || isempty(tail)
    tail                        = 0;
elseif isnumeric(tail) && isscalar(tail) && ismember(tail,[-1 0 1])
    % OK, grandfathered
else
    [ ~, tail ]                 = internal.stats.getParamVal(tail,{'left','both','right'},'TAIL');
    tail                        = tail - 2;
end

samplesize                      = nx;
xmean                           = mx;
ser                             = sx;
df                              = max(samplesize - 1,0);
tval                            = ( xmean - m ) ./ ser;
if nargout > 3
    stats                       = struct('tstat', tval, 'df', cast(df,class(tval)), 'sd', sdpop);
    if isscalar(df) && ~isscalar(tval)
        stats.df                = repmat(stats.df,size(tval));
    end
end

% Compute the correct p-value for the test, and confidence intervals
% if requested.
if tail == 0 % two-tailed test
    p                           = 2 * tcdf(-abs(tval), df);
    if nargout > 2
        crit                    = tinv((1 - alpha / 2), df) .* ser;
        ci                      = cat(dim, xmean - crit, xmean + crit);
    end
elseif tail == 1 % right one-tailed test
    p                           = tcdf(-tval, df);
    if nargout > 2
        crit                    = tinv(1 - alpha, df) .* ser;
        ci                      = cat(dim, xmean - crit, Inf(size(p)));
    end
elseif tail == -1 % left one-tailed test
    p                           = tcdf(tval, df);
    if nargout > 2
        crit                    = tinv(1 - alpha, df) .* ser;
        ci                      = cat(dim, -Inf(size(p)), xmean + crit);
    end
end
% Determine if the actual significance exceeds the desired significance
h                               = cast(p <= alpha, class(p));
h(isnan(p))                     = NaN; % p==NaN => neither <= alpha nor > alpha

return % local_ttesthat

%------------------------------------------------------------------------
% p = local_ttest2hat( mx, sx, nx, my, sy, ny )
% identical to ttest2, works on summary stats (rather than raw data)
%------------------------------------------------------------------------
function p = local_ttest2hat( mx, sx, nx, my, sy, ny )

s2x                             = nx .* ( sx .^ 2 );
s2y                             = ny .* ( sy .^ 2 );
difference                      = mx - my;
s2xbar                          = s2x ./ nx;
s2ybar                          = s2y ./ ny;
dfe                             = ( s2xbar + s2ybar ) .^2 ./ ( s2xbar.^2 ./ (nx-1) + s2ybar.^2 ./ (ny-1));
se                              = sqrt(s2xbar + s2ybar);
ratio = difference ./ se;
if se == 0
    dfe                         = 1; 
end
p                               = 2 * tcdf(-abs(ratio),dfe);

return

%------------------------------------------------------------------------
% [ idx, vals ] = local_find_local_extrema( x, mode )
%  detect all local extrema in a matrix
%------------------------------------------------------------------------
function [ idx, vals ] = local_find_local_extrema( x, mode )

if nargin < 2 || isempty( mode )
    mode                        = '';
end
x                               = x( : );
[ m, n ]                        = size( x );
d2                              = diff( sign( diff( x ) ) );
switch mode
    case 'ext'
        [ row, col ]            = find( abs( d2 ) > 1 );
    case 'min'
        [ row, col ]            = find( d2 > 1 );
    otherwise
        [ row, col ]            = find( d2 < -1 );
end
row                             = row + 1;
if n == 1
    idx                         = row;
else
    idx                         = [ row col ];
end
vals                            = x( row + ( col - 1 ) * m );

return % local_find_local_extrema

% EOF