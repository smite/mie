% Mie3Dplot  a 3D polar plot of Mie-scattered irradiance
%   R = Mie3Dplot(x, m [, N, [, s]])
%   Input:  x is the size parameter of the sphere,
%     m is the refractive index of the scattering sphere
%     relative to medium, N the required angular precision and
%     s the scalar surface conductance parameter
%   Output: Array R, size(R) = [l(x), l(m), l(theta) = N+1, l(phi) = N/2+1]

% Calls MieStokes.m
% Ville Bergholm 2002



function R = Mie3Dplot(x, m, N, s)

if nargin < 4
  s = 0;  % no surface conductance as a default
  if nargin < 3
    N = 100;
  end
end

theta = (0:N)'*pi/N;     % from zero to pi 
phi = (0:2:N)'*pi/(2*N); % col. vector, from zero to pi/2
I = [ones(length(phi),1), cos(2*phi), sin(2*phi), zeros(length(phi),1)];
O = MieStokes(x, m, theta, I, s);

% just irradiance, arbitrary growing function applied
R = log10(O(:,:,:,:,1) + 1);

%R = R - min(min(R));

stheta = sin(theta)*ones(1, length(phi));
ctheta = cos(theta)*ones(1, length(phi));

sphi = ones(length(theta),1)*sin(phi');
cphi = ones(length(theta),1)*cos(phi');

lx = length(x);
lm = length(m);

% loop here
for i = 1:lx
   for j = 1:lm

      temp = squeeze(R(i,j,:,:)).*stheta;
      z = squeeze(R(i,j,:,:)).*ctheta;
      x = temp.*cphi;
      y = temp.*sphi;

      figure((i-1)*lm+j);
      surf(x, y, z, squeeze(R(i,j,:,:)));
      ots1 = 'Angular distribution of scattered irradiance, incoming wave x-polarized. x=';
      ots2 = ', m=';
      xlabel('x');
      ylabel('y');
      zlabel('z');
      view([-60, 20]) % useful viewing angle
      shading interp  % gouraud shading of colours
      lighting phong  % nice
      axis equal      % equal scales on axes
      axis on         % draw axes
      light('Position',[-1 -2 0.6]);
      %title(strcat(ots1, num2str(x,3), ots2, num2str(m,3)));
   end
end
