% CALC_BIAS_PT          calculate an analytical estimate of MI bias from a joint probability table.
%
% call                  BIAS = CALC_BIAS_PT( MAT, MODE )
%
% gets                  MAT         matrix - rows are "responses" (activity)
%                                   and columns are "stimuli" (behaviors)
%                       MODE        0 or 'none'     nothing
%                                   1 or 'simple'   underestimate
%                                   2 or 'bayes'    supposed to be optimal
%                                   3 or 'full'     overestimate
%
% return                BIAS        in bits
%
% calls                 BAYESCOUNT
%
% calles by             DECODEPOP, INFO_CALCULATOR, MUTINF_NEW

% REFERENCES
% 1. Panzeri and Treves 1996, Network
% 2. Panzeri et al., 2007, JNP review
% 3. BAYESCOUNT was downloaded on 06-aug-07 from 
%   http://stefano.panzeri.googlepages.com/informationbiascorrections

% 11-aug-07 ES

function bias = calc_bias_PT( mat, mode )

N = sum( sum( mat ) );
S = size( mat, 2 );
switch mode
    case { 0, 'none' }
        bias = 0;
        return
    case { 1, 'simple' } % may underestimate
        sumR_s = sum( sum( mat ~= 0 ) );
        R = sum( sum( mat, 2 ) ~=0 );
    case { 2, 'bayes' } % supposed to be optimal
        R = bayescount( sum( sum( mat ) ), sum( mat, 2 ) / sum( sum( mat ) ) );
        nz = sum( mat );
        R_s = zeros( 1, S );
        for s = 1 : S
            if nz( s )
                R_s( s ) = bayescount( sum( mat( :, s ) ), mat( :, s ) / sum( mat( :, s ) ) );
            end
        end
        sumR_s = sum( R_s );
    case { 3, 'full' } % the largest possible bias - may overestimate
        sumR_s = prod( size( mat ) );
        R = size( mat, 1 );
    otherwise
        error( 'unsupported mode' )
end
bias = ( sumR_s - R - S + 1 ) / ( 2 * N * log( 2 ) );

% the previousy used version was an inline, equivalent to 'simple' mode:
% calc_biasPT = inline( ...
%     '( sum(sum(m~=0)) - sum(sum(m,2)~=0) - size(m,2) + 1 ) / ( 2 * sum(sum(m)) * log( 2 ) )'...
%     , 'm' );

return
