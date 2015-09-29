% MiePlot_xth(m [, N])
% Makes surface plots of Mie scattered irradiances
% as a fuction of size parameter x and scattering angle theta.
% m = k_p / k_m is the ratio refractive indices.
% N is the number of grid points in both directions.

% Ville Bergholm 2002

function MiePlot_xth(m, N)

if nargin < 2
   N = 80;
end

step = (15-1)/N;
x=1:step:15;
step = 180/N;
theta=0:step:180;
th=theta*pi/180;

[i1,i2,P] = MieSIrr(x,m,th);

R1 = squeeze(log10(i1));
R2 = squeeze(log10(i2));
R3 = squeeze(log10(i1+i2));

figure(1);
surf(theta, x, R1);

zlabel('log_{10}(I_s / I_i)');
ylabel('Size parameter x');
xlabel('Scattering angle \theta');
title(strcat('Perpendicular scattered irradiance, m=', num2str(m,3)));
view([22, 28]) % useful viewing angle
shading interp  % gouraud shading of colours
lighting phong  % nice
% axis equal      % equal scales on axes
axis on         % draw axes
%light('Position',[-1 -2 0.6]);
light('Position',[130 -6 3]);

figure(2);
surf(theta, x, R2);

zlabel('log_{10}(I_s / I_i)');
ylabel('Size parameter x');
xlabel('Scattering angle \theta');
title(strcat('Parallel scattered irradiance, m=', num2str(m,3)));
view([22, 28]) % useful viewing angle
shading interp  % gouraud shading of colours
lighting phong  % nice
% axis equal      % equal scales on axes
axis on         % draw axes
%light('Position',[-1 -2 0.6]);
light('Position',[130 -6 3]);


%return;

figure(3);
surf(theta, x, R3);
zlabel('log_{10}(I_s / I_i)');
ylabel('Size parameter x');
xlabel('Scattering angle \theta');
title(strcat('Total scattered irradiance, m=', num2str(m,3)));
view([22, 28]) % useful viewing angle
shading interp  % gouraud shading of colours
lighting phong  % nice
% axis equal      % equal scales on axes
axis on         % draw axes
%light('Position',[-1 -2 0.6]);
light('Position',[130 -6 3]);

%return;

figure(4);
surf(theta, x, P);
ylabel('Size parameter x');
xlabel('Scattering angle \theta');
title(strcat('Polarization, m=', num2str(m,3)));
view([22, 28]) % useful viewing angle
shading interp  % gouraud shading of colours
lighting phong  % nice
% axis equal      % equal scales on axes
axis on         % draw axes
%light('Position',[-1 -2 0.6]);
light('Position',[130 6 5]);
