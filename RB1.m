% RB1   Riccati-Bessel function of the first kind.
%       psi = RB1(nmax, rho [, acc]), where
%       psi(n, m) = psi_n(rho(m)) = rho(m) * j_n(rho(m))
%       n = 2..nmax, m = 1..length(rho)

% Uses the method from W. J. Lentz, "Generating Bessel functions in Mie scattering calculations using continued fractions," Appl. Opt. 15, 668- (1976)
% Ville Bergholm 2001-2008

function psi = RB1(nmax, rho, acc)

if nargin < 3
  acc = 0.000001; % desired accuracy
end

% rho into a row vector
rho = rho(:).';

% Lentz initialization: we shall calculate the "exact"
% ratio j_{nmax-1}/j_nmax using a continuous fraction

n = nmax;    % a little shorthand...
ratio = 2*(n+0.5)./rho; % a1 == b1
a(1,:) = ratio;
b(1,:) = ratio;
c(1,:) = Inf*ones(size(rho)); % trick, 1/c(1,:) = 0
s=1;
% ratio will be  J_{n-0.5} / J_{n+0.5} = \psi_{n-1} / \psi_n

while 1
  a(s+1,:) = (-1)^s*2*(n+0.5+s)./rho;
  b(s+1,:) = 1./b(s,:) + a(s+1,:);
  c(s+1,:) = 1./c(s,:) + a(s+1,:);
  s=s+1;
  if max(abs(b(s,:)-c(s,:))) < acc
    break;
  else
    ratio = ratio.*b(s,:)./c(s,:);
  end;
end;

% now we'll use the ratio to intialize downward recursion:

psi = zeros(nmax, length(rho));
psi(nmax,:) = 1e-10;
psi(nmax-1,:) = psi(nmax,:).*ratio;

%  psi(1,:) = sin(rho)./rho - cos(rho);
%  psi(2,:) = 3*psi(1,:)./rho - sin(rho);

% upward recursion is unstable, we'll use downward recursion
for n=nmax-2:-1:1
  psi(n,:) = (2*n+3)*psi(n+1,:)./rho - psi(n+2,:);
end

% We know that psi_0 == sin(rho), use this to scale all psi values
fix = sin(rho) ./ (3*psi(1,:)./rho - psi(2,:));

psi = psi.*(ones(nmax,1) * fix);
