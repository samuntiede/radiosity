function CDF = Scaled_BetaCDF(y, a, b, p, q)
%   Scaled_BetaCDF - General (4-paramter) Beta cumulative density function.
%
%   Y = Scaled_BetaCDF(y, a, b, p, q) returns the beta probability density
%   function with parameters a and b at the values in y, where the minimum
%   and maximum limits of the distribution are given by p and q,
%   respectively
%
% Usage:
% CDF = Scaled_BetaCDF(y, a, b);
% CDF = Scaled_BetaCDF(y, a, b, p);
% CDF = Scaled_BetaCDF(y, a, b, p, q);
% 
%
% See also: Scaled_BetaPDF, Scaled_BetaStat
%
%
% Author: Greig A. Paterson
%

if nargin < 3
    error('Scaled_BetaCDF:Input', 'At least 3 inputs are required.')
end

% Default to standard beta distribution
if nargin < 4
    p = 0;
end

if nargin < 5
    q = 1;
end


if a <= 0 || b <= 0
    error('Scaled_BetaCDF:Shape', 'Shape paramers must be greater than zero.')
end

if q <= p
    error('Scaled_BetaCDF:Range', 'Invalid range, q must be greater than p.')
end

% intialize the CDF to zeros
CDF = NaN(size(y));

% Set probabilities to zero outside of range
CDF(y < p) = 0;
CDF(y > q) = 1;

% Create dummy variable to handle in range values only
k = find(isnan(CDF));
x = y(k);


CDF(k) = betainc((x-p)./(q-p),a,b);


            
