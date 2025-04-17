function [M, V, Med, Mo, Sk, exK] = Scaled_BetaStat(a, b, p, q)
%   Scaled_BetaStat - General (4-paramter) Beta statistics.
%
%   [M, V, Med, Mo, Sk, K] = Scaled_BetaStat(y, a, b, p, q) returns the 
%   parametric statisitcs for a general beta distribution with parameters 
%   a and b and the minimum and maximum limits of the distribution given 
%   by p and q, respectively
%   
%   M - Mean
%   V - Variance
%   Med - Median
%   Mo - Mode
%   Sk - Skewness
%   exK = Excess Kurtosis (normal distribution has exK == 0)
%
%
% Usage:
% M = Scaled_BetaStat(y, a, b);
% M = Scaled_BetaStat(y, a, b, p);
% M = Scaled_BetaStat(y, a, b, p, q);
% [M, V, Med, Mo, Sk, exK] = Scaled_BetaStat(a, b, p, q);
% 
%
% See also: Scaled_BetaCDF, Scaled_BetaPDF
%
%
% Author: Greig A. Paterson
%


if nargin < 2
    error('Scaled_BetaStat:Input', 'At least 2 inputs are required.')
end

% Default to standard beta distribution
if nargin < 3
    p = 0;
end

if nargin < 4
    q = 0;
end

if a <= 0 || b <= 0
    error('Scaled_BetaCDF:Shape', 'Shape paramers must be greater than zero.')
end

if q <= p
    error('Scaled_BetaCDF:Range', 'Invalid range, q must be greater than p.')
end


% The mean
M = (a.*q + b.*p) ./ (a + b) ;

% The varaince 
V = ( a.*b.*(q - p).^2 ) ./ ( (a+b).^2 .* (a + b + 1) ) ;

% The Median
Med = p + (q - p) .* (betaincinv(0.5,a,b));

% The Mode
if a > 1 && b > 1 
    Mo = ( q.*(a-1) + p.*(b-1) ) ./ (a+b-2);
else
    Mo=NaN;
end

% The skewness
Sk = ( 2.*(b - a).*sqrt(a + b + 1) ) ./ ( (a + b + 2).*sqrt(a.*b) );

% The kurtosis
exK = ( 6.*( (a - b).^2 .* (a + b + 1) - a.*b.*(a + b + 2) ) ) ./ (a.*b.*(a + b + 2).*(a + b + 3));


