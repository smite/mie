function [S1, S2] = MieSMatrix(x, m, theta, s)
% MieSMatrix  Mie scattering matrix elements S1 and S2.
%       [S1, S2] = MieSMatrix(x, m, theta [, s])
%       Input:
%         (v) x = k_m*a is the dimensionless size parameter,
%         (v) m = k_p/k_m is the refractive index of the particle relative to the medium,
%         (v) theta is the scattering angle, and
%         (s) s = sigma_s*c_m*mu_m is the optional surface conductivity parameter.
%       Output:
%         Arrays S1 and S2, size(S1) = [length(x), length(m), length(theta)]

% Ville Bergholm 2001-2008

if nargin < 4
  s = 0;
end

nmax = MieWn(x);
disp('MiePT begins');
tic
[p, t] = MiePT(nmax, theta); % rows: n, columns: theta
toc
disp('MieAB begins');
tic
[a, b] = MieAB(nmax, x, m, s); % n, x, m
toc

% TODO kaikille x:ille oma nmax, loput termit nollaksi?
disp('NaN check begins');
tic
% modified Barnett's check for NaN:s
% Check for invalid (NaN) results due to too many terms in
% relatively small particles.
while 1   % infinite loop
   test = find(any(any(isnan([a;b]), 3),1)); % indices of faulty x's
   if isempty(test)
      break; % break out of the while loop
   end
   disp('MieSc: NaN found');
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
toc

disp('precalculation begins');
tic
n=(1:nmax).'; % column vector

% precalculated numerical factor
%E = 0.5*((2*n+1)./(n.*(n+1)))*ones(1,length(x));
%a = E.*a;
%b = E.*b;
%Sp = ((a+b).')*(p+t);
%Sm = ((a-b).')*(p-t);

E = 0.5*((2*n+1)./(n.*(n+1)))*ones(1,length(theta));
cp = (a+b);
cm = (a-b);
qp = E.*(p+t);
qm = E.*(p-t);
toc
disp('contraction starts');
tic
% if Matlab had a proper contraction operator, all this wouldn't be needed
%lx = length(x);
%lm = length(m);
%ltheta = length(theta);
%Sp = zeros(lx,lm,ltheta);
%Sm = Sp;
%for j=1:lx
%   for k=1:lm
%      for l=1:ltheta
%         for i=1:nmax
%            Sp(j,k,l) = Sp(j,k,l) + cp(i,j,k)*qp(i,l);
%            Sm(j,k,l) = Sm(j,k,l) + cm(i,j,k)*qm(i,l);
%         end
%      end
%   end
%end

% mex-file "contract.dll" replaces the preceding loops (about 100 times faster!)
Sp = contract(cp,qp,1,1);
Sm = contract(cm,qm,1,1);

% up to here, Sp(j,k,l) = cp(*,j,k)qp(*,l)
toc

S1 = Sp+Sm;
S2 = Sp-Sm;
