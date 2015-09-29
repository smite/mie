% MieLobe  Lobe clusters in Mie-scattered light
%    [i1, i2, P] = MieLobe(x, m, thc, thr, N [,s]) returns the
%   Input: Size parameter x = k_m*a, refractive index of particle
%     relative to medium m = k_p / k_m, theta center thc, theta radius thr,
%     number of samples N, scalar surface conductance parameter s
%   Output: Arrays i1, i2 and P, size(i1) = [length(x), length(m)].
%     These contain the relative scattered irradiance of
%     incident light perpendicular (i1) and parallel (i2) to the
%     scattering plane.
%     P gives the polarization of the scattered light.

%   Ville Bergholm 2002-2008

function [i1, i2, P] = MieLobe(x, m, thc, thr, N, s)

if nargin < 6
  s = 0;  % no surface conductance as a default
end

%theta = thc + (-1 + 2*(1:N)/(N+1))*thr;
theta = thc + (-1 + 2*(0:N)/N)*thr;

[S1, S2] = MieSMatrix(x, m, theta, s);

a1 = zeros(length(x), length(m));
a2 = a1;

coef = 1/(2*N+1);
% average over theta values

for i = 1:(2*N+1) %loop over theta values
   a1(:,:) = a1(:,:) + S1(:,:,i)*coef;
   a2(:,:) = a2(:,:) + S2(:,:,i)*coef;
end

i1 = a1.*conj(a1);
i2 = a2.*conj(a2);
P = (i1-i2)./(i1+i2);

%axis auto;
%xlabel('Size parameter x');
%ylabel('lg(I_s/I_i)');
%plot(x, i1, 'b-',x, i2, 'r-');
%legend('i1, perpendicular', 'i2, parallel', 0);
