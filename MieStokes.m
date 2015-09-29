function O = MieStokes(x, m, theta, I, s)
% MieStokes  Gives the Stokes' vector O of Mie-scattered light.
%    O = MieStokes(x, m, theta, I [, s])
%    Input:
%      (v) x = k_m*a is the dimensionless size parameter,
%      (v) m = k_p/k_m is the refractive index of the particle relative to the medium,
%      (v) theta is the scattering angle,
%      (m) I = (I, Q, U, V;...) is a matrix with Stokes' vectors of the incident light, and
%      (s) s = sigma_s*c_m*mu_m is the optional surface conductivity parameter.
%    Output:
%      Array O, size(O) = [lenght(x), length(m), length(theta), size(I,1), 4]

% Ville Bergholm 2001-2008

if nargin < 5
   s = 0;
end

[S1, S2] = MieSMatrix(x, m, theta, s);

S11 = 0.5*(S1.*conj(S1) + S2.*conj(S2));
S12 = 0.5*(S2.*conj(S2) - S1.*conj(S1));
S33 = 0.5*(S1.*conj(S2) + S2.*conj(S1));
S34 = 0.5i*(S1.*conj(S2) - S2.*conj(S1));
% size(Sxy) = [length(x), length(m), length(theta)]

% calculate the output stokes' vectors
n = size(I,1);
O = zeros(length(x), length(m), length(theta), n, 4);
for i = 1:n
   O(:,:,:,i,1) = I(i,1)*S11 + I(i,2)*S12;
   O(:,:,:,i,2) = I(i,1)*S12 + I(i,2)*S11;
   O(:,:,:,i,3) = I(i,3)*S33 + I(i,4)*S34;
   O(:,:,:,i,4) = -I(i,3)*S34 + I(i,4)*S33;
end

%O=squeeze(O); % remove singleton dimensions
