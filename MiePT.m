function [p, t] = MiePT(nmax, theta)
% MiePT  Angle-dependent Mie scattering functions pi and tau.
%        \pi_n(\theta) := P_n^1 (\cos \theta) / (\sin \theta),
%        \tau_n(\theta) := d P_n^1 (\cos \theta) / (d \theta),
%        where P_n^m is the associated Legendre polynomial.
%
%        [pi, tau] = MiePT(nmax, theta)
%        Input:
%          nmax is the highest required degree of pi and tau,
%          theta is a vector of scattering angles.
%        Output:
%          arrays p and t, size(p) = [nmax, length(theta)],
%          p[n, i] = \pi_n(theta(i)).

% Uses recursion formulas from W.J. Wiscombe, "Improved Mie scattering algorithms", Applied Optics 19, 1505- (1980).
% Ville Bergholm 2001-2008

theta = theta(:).'; % theta is a row vector
temp = cos(theta);

p(1,:) = ones(1, length(theta));
t(1,:) = temp;
p(2,:) = 3*temp;
for n = 2:nmax-1
  x = temp.*p(n,:); % temp variables, called s and t in Wiscombe's paper
  y = x - p(n-1,:);
  t(n,:) = n*y - p(n-1,:);
  p(n+1,:) = x + ((n+1)/n)*y;

  % "traditional" recursion formulas:
  % p(n,:) = ((2*n-1)*temp.*p(n-1,:) - n*p(n-2,:))/(n-1);
  % t(n,:) = n*temp.*p(n,:) - (n+1)*p(n-1,:);
end

t(nmax,:) = nmax*temp.*p(n,:) - (nmax+1)*p(nmax-1,:);
