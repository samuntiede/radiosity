function PDF = Scaled_BetaPDF(y, a, b, p, q)
%   Scaled_BetaPDF - General (4-paramter) Beta probability density function.
%
%   Y = Scaled_BetaPDF(y, a, b, p, q) returns the beta probability density
%   function with parameters a and b at the values in y, where the minimum
%   and maximum limits of the distribution are given by p and q,
%   respectively
%
% Usage:
% PDF = Scaled_BetaPDF(y, a, b);
% PDF = Scaled_BetaPDF(y, a, b, p);
% PDF = Scaled_BetaPDF(y, a, b, p, q);
% 
%
% See also: Scaled_BetaCDF, Scaled_BetaStat
%
%
% Author: Greig A. Paterson
%

if nargin < 3
    error('Scaled_BetaPDF:Input', 'At least 3 inputs are required.')
end

% Default to standard beta distribution
if nargin < 4
    p = 0;
end

if nargin < 5
    q = 1;
end

if a <= 0 || b <= 0
    error('Scaled_BetaPDF:Shape', 'Shape paramers must be greater than zero.')
end

if q <= p
    error('Scaled_BetaPDF:Range', 'Invalid range, q must be greater than p.')
end

PDF = ( (y-p).^(a-1) .* (q - y).^(b-1) ) ./ ( (q - p).^(a+b-1) .* beta(a,b) );

% Ensure probabilities to zero outside of range
PDF(y < p) = 0;
PDF(y > q) = 0;