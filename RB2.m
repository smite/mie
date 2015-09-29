% RB2   Riccati-Bessel function of the second kind.
%       zeta = RB2(nmax, rho), where
%       zeta(n, m) = zeta_n(rho(m)) = rho(m) * h_n^{(1)}(rho(m))
%       n = 2..nmax, m = 1..length(rho)

% Ville Bergholm 2001-2008

function zeta = RB2(nmax, rho)

zeta = zeros(nmax, length(rho));
zeta(1,:) = -cos(rho)./rho - sin(rho);
zeta(2,:) = 3*zeta(1,:)./rho + cos(rho);

% upward recursion is stable
for n=3:nmax
  zeta(n,:) = (2*n-1)*zeta(n-1,:)./rho - zeta(n-2,:);
end
