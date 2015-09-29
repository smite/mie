function n = n_Air(t, p)
% n_Air  Index of refraction for dry air at 630 nm wavelength (He-Ne laser)
%        n = n_Air(T [, p])
%        T is temperature in degrees Celsius.
%        p is pressure in Pascals. 101325 Pa is the default.

% Formula from CRC
% Ville Bergholm 2002-2008

if nargin < 2
   p = 101325;
end

% lambda = 630 nm (He-Ne)
% r = (n-1)*1e8 
r = 27656;

% r = 27728; % lambda = 580 nm

% magical formula from CRC
n = 1 + 1e-8*r*p.*(1 + p*(61.3 - t)*1e-10)./(96095.4*(1 + 0.003661*t));
