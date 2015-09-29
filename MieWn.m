function nmax = MieWn(x)
% MieWn  Wiscombe approximation for the number of terms in Mie series.
%        nmax = MieWn(x), where x = k_m*a is the dimensionless size parameter.
%        Only works for real x.

% Uses the heuristic from W.J. Wiscombe, "Improved Mie scattering algorithms", Applied Optics 19, 1505- (1980).
% Ville Bergholm 2002-2008

x = x(:);

if (any(imag(x)))
  error('x must be real.')
end

temp = max(abs(x));
if temp < 0.02
  error('x is too small.')
elseif temp <= 8
  nmax = ceil(temp+4*(temp^(1/3))+1)
elseif temp < 4200
  nmax = ceil(temp+4.05*(temp^(1/3))+2)
elseif temp <= 20000
  nmax = ceil(temp+4*(temp^(1/3))+2)
else
  error('x is too large.')
end
