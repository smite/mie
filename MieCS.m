% MieCS  Mie scattering and extinction cross sections.
%   [Cext, Csca] = MieCS(a, k_m, m [, s])
%   Input: Scattering particle radius a,
%     wave vector of scattering light in surrounding medium k_m,
%     refractive index of particle relative to medium m = k_p / k_m
%     and the optional surface conductivity parameter s = sigma_s * c_m * mu_m.
%   Output: Arrays Cext and Csca, size(Cext) = [length(x), length(m)]

% Ville Bergholm 2001-2008

function [Cext, Csca] = MieCS(r, k, m, s)

if nargin < 4
  % particle surface conductivity parameter
  s = 0;
end

r = r(:);
k = k(:);
m = m(:);
x = k.*r; % no tricks here, keep it simple... should make this more general 

% Wiscombe approximation for # of terms in series
% only works for real x
nmax = MieWn(x);

[a,b] = MieAB(nmax, x, m, s);
%size(a)

% TODO also here have an nmax for each x separately (in MieAB?)
% Are there NaN:s among the results?
while 1 %infinite loop
   test = find(any(any(isnan([a;b]), 3),1)); % indices of faulty x's
   if isempty(test)
      break; % break out of the while loop
   end
   disp('MieCS: NaN found');
   a(:, test, :) = 0;
   b(:, test, :) = 0; % remove the NaN:s
   nmax2 = MieWn(x(test)); % take a new (smaller) nmax
   if nmax2 > nmax
      error('This should never happen');
   end
   [A,B] = MieAB(nmax2, x(test), m); %compute new A, B for the invalid x values
   a(1:nmax2, test, :) = A;
   b(1:nmax2, test, :) = B; % splice the results into a, b
end

n=1:nmax;
N = 2*n+1;
%size(N)

ex = real(a+b);
sc = abs(a).^2 + abs(b).^2;
%size(sc)

% shiftdim gets rid of leading singleton dimension (from N)
Cext = shiftdim(contract(N, ex, 2, 1), 1);
Csca = shiftdim(contract(N, sc, 2, 1), 1);
%size(Csca) % = [l(x), l(m)]

% either r or k can be vectors
temp = (2*pi)./(k.^2);
temp = temp.*ones(length(r),1);
temp = temp*ones(1,length(m));
%size(temp)

Cext = temp.*Cext;
Csca = temp.*Csca;
%size(Cext)
