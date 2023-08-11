function [AB_pt, AB_sim] = SurrogateAcceptanceBands(c0, c_jit, alpha)
% Generates pointwise and simultaneous acceptance bands from an uncorrected
% CCG (input c0) and jittered or shuffled surrogate CCGs (input cm_all).
% See section 3.1 of https://arxiv.org/pdf/1111.4296.pdf
%
% The number of shuffle/jitter repetitions, size(cm_jit,2), should be
% divisible by alpha/2 for indexing purposes. If it is not, then alpha
% will be adjusted to a value displayed in a warning message.
%
% INPUTS
% - c0: [time lags x 1] uncorrected CCG
% - c_jit: [time lags x jitter repetitions] matrix where each column is a
%       CCG generated from a different shuffle/jitter surrogate except for
%       the first, which should be the original uncorrected/unjittered CCG.
% - alpha: significance level. E.g., 0.05 for 95% acceptance bands.
%
% OUTPUTS
% - AB_pt: [time lags x 2] pointwise acceptance band. Col 1 is lower, col 2
%       is upper.
%       * Can be used for hypothesis testing at one time lag (tau) chosen a
%         priori, probably tau = 0 to test whether two neurons are firing
%         synchronously.
% - AB_sim: [time lags x 2] simultaneous acceptance band where Col 1 is
%       lower and col 2 is upper.
%       * This one is corrected for multiple comparisons using a 
%         proprietory blend of the finest z-scoring(?) and ordinal stats.

M = size(c_jit, 2); % number of surrogate CCGs
c = [c0 c_jit];
c_ord = sort(c, 2); % for each lag t, sort all values in ascending order

% Calculate positions in the ordered arrays to get as lower/upper bounds.
% +1 for both because Amarasingham et al. index from 0 in their
% mathematical notation and MATLAB indexes from 1.
try
    a_ind = M*(alpha/2) + 1;
    b_ind = M*(1-alpha/2) + 1;
catch
    a_ind = floor(M*(alpha/2)) + 1;
    b_ind = ceil(M*(1-alpha/2)) + 1;
    new_alpha = ((a_ind-1)/M)*2;
    new_percent = 1 - new_alpha;
    warning('Alpha = %d and number of shuffle/jitter iterations M = %d are incompatible for array indexing!\nUsing alpha = %d instead to generate %d percent acceptance bands.', alpha,M,new_alpha,new_percent);
end

%% pointwise acceptance bands
AB_pt = [c_ord(:,a_ind), c_ord(:,b_ind)];

%% simultaneous acceptance bands
% this part exactly mirrors a block of equations in the paper
v = mean(c_ord(:,2:end-1), 2); % robust to one extreme outlier in each tail, e.g. if there is a significant peak in the CCG
s = std(c_ord(:,2:end-1), 0, 2);
c_star = (c - v) ./ s; % equivalent to z-scoring c, but ignoring the single largest and smallest values in the mean and std calculations
c_star_max_ord = sort(max(c_star,[],1));
c_star_min_ord = sort(min(c_star,[],1));

a_star = (c_star_min_ord(a_ind) * s) + v;
b_star = (c_star_max_ord(b_ind) * s) + v;

AB_sim = [a_star b_star];

end

