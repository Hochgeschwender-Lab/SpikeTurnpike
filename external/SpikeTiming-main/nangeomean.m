function m = nangeomean(x,dim)
%NANGEOMEAN Geometric mean.
% 
% for more details, see geomean.m
% two differences:  (1) NaNs are excluded from the averaging
%                   (2) zeros are replaced by machine precision ("eps")

% 18-mar-13 ES

% revisions
% 10-dec-19 added eps to prevent m=0 when there is a single 0 in x

if any(x(:) < 0)
    error(message('stats:geomean:BadData'))
end

if nargin < 2 || isempty(dim)
    % Figure out which dimension sum will work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

n = size(x,dim);
% Prevent divideByZero warnings for empties, but still return a NaN result.
if n == 0, n = NaN; end

% Take the n-th root of the product of elements of X, along dimension DIM.
% if nargin < 2
%     m = exp(nanmean(log(x+eps)));
% else
%     m = exp(nanmean(log(x+eps),dim));
% end
if nargin < 2
    m = exp(nansum(log(x+eps))./n);
else
    m = exp(nansum(log(x+eps),dim)./n);
end