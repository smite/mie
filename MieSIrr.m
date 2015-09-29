function [i1, i2, P] = MieSIrr(x, m, theta, s)
% MieSIrr  Mie-scattered irradiance and polarization.
%    [i1, i2, P] = MieSIrr(x, m, theta [, s])
%    Input:
%      (v) x = k_m*a is the dimensionless size parameter,
%      (v) m = k_p/k_m is the refractive index of the particle relative to the medium,
%      (v) theta is the scattering angle, and
%      (s) s = sigma_s*c_m*mu_m is the optional surface conductivity parameter.
%    Output:
%      Arrays i1, i2 and P, size(i1) = [length(x), length(m), length(theta)].
%      These contain the relative scattered irradiance of
%      incident light perpendicular (i1) and parallel (i2) to the
%      scattering plane as a function of the scattering angle.
%      P gives the polarization of the scattered light.

% Ville Bergholm 2001-2008

if nargin < 4
  s = 0;  % no surface conductance as a default
end

[S1, S2] = MieSMatrix(x, m, theta, s);

i1 = S1.*conj(S1);
i2 = S2.*conj(S2);
P = (i1-i2)./(i1+i2);

%axis auto;
%xlabel('Theta (rad)');
%ylabel('lg(I_s/I_i)');
%plot(theta, squeeze(i1), 'b-', theta, squeeze(i2), 'r-');
%legend('i1, perpendicular', 'i2, parallel', 0);
