function [a,b] = MieAB(nmax, x, m, s)
% MieAB  Mie Scattering coefficients a_i and b_i.
%        [a, b] = MieAB(nmax, x, m, [, s])
%        Input:
%          nmax is the highest required order of a and b,
%          x = k_m*a is the dimensionless size parameter,
%          m = k_p/k_m is the refractive index of the particle relative to the medium,
%          and s = sigma_s*c_m*mu_m is the optional surface conductivity parameter.
%        Output:
%          Arrays a and b, size(a) = [nmax, length(x), length(m)],
%          a(n, i, j) = a_n(x(i), m(j))

%   Ville Bergholm 2002-2008

if nargin < 4
  s = 0; % default: no surface conductance
end

% open up m and x into row vectors
m = m(:).';
x = x(:).';
mx = kron(m,x); % long vector

% Riccati-Bessel functions
psi = RB1(nmax, x);     % x is a row vector, size(psi) = [nmax, length(x)] 
psi_m = RB1(nmax, mx);% size = [nmax, length(mx)]
zeta = RB2(nmax, x);
xi = psi + 1i * zeta;  % R-B function corresponding to 1st sph. Hankel


% temp vectors for calculating derivatives
Tpsi = [sin(x) ; psi(1:(nmax-1),:)];     % psi_0 = sin(x)
Tpsi_m = [sin(mx) ; psi_m(1:(nmax-1),:)]; %longer
Tzeta = [-cos(x) ; zeta(1:(nmax-1),:)];  % zeta_0 = -cos(x)

% temp matrices to facilitate computation of derivatives
N = ((1:nmax).')*ones(1,length(x));
Nm = ((1:nmax).')*ones(1,length(mx));

if length(mx) > 1
   mx2 = ones(nmax,1)*mx;
   if length(x) > 1
      x2 = ones(nmax,1)*x;
   else
      x2 = x;
   end
else
   mx2 = mx;
   x2 = x;
end

% derivatives are calculated using recursion formulae
Dpsi = Tpsi-N.*psi./x2;
Dpsi_m = Tpsi_m-Nm.*psi_m./(mx2); %
Dzeta = Tzeta-N.*zeta./x2;
Dxi = Dpsi + 1i * Dzeta;

temp = ones(1,length(m)); % to make use of the kronecker products
C = 1i*s.*temp; %s has to be scalar

% finally, Mie coefficients
%a = (m.*psi_m.*Dpsi - psi.*Dpsi_m + C.*Dpsi.*Dpsi_m) ./ (m.*psi_m.*Dxi - xi.*Dpsi_m + C.*Dpsi_m.*Dxi);
%b = (psi_m.*Dpsi - m.*psi.*Dpsi_m + C.*psi.*psi_m) ./ (psi_m.*Dxi - m.*xi.*Dpsi_m + C.*psi_m.*xi);

% works as long as in kron(a,b) a is vector, b is matrix
a = (psi_m.*kron(m,Dpsi) - kron(temp,psi).*Dpsi_m + kron(C,Dpsi).*Dpsi_m) ./ (psi_m.*kron(m,Dxi) - kron(temp,xi).*Dpsi_m + kron(C,Dxi).*Dpsi_m);
b = (psi_m.*kron(temp,Dpsi) - kron(m,psi).*Dpsi_m + kron(C,psi).*psi_m) ./ (psi_m.*kron(temp,Dxi) - kron(m,xi).*Dpsi_m + kron(C,xi).*psi_m);
shape = [nmax, length(x), length(m)]; % reshape the kron-deformed array
a = reshape(a, shape);
b = reshape(b, shape);
