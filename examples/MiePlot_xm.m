% MiePlot_xm(theta, [, N])
% Makes surface plots of Mie scattered irradiances
% as a fuction of parameters x and m.
% Theta is the scattering angle in degrees.
% N is the number of grid points in both directions.

% Ville Bergholm 2002

function MiePlot_xm(theta, N)

if nargin < 2
   N = 80;
end

step = (15-1)/N;
x=1:step:15;
% remember a little offset to avoid making m exactly zero
step = (3-0.09)/N;
m=0.09:step:3;
th=theta*pi/180;

[i1,i2,P] = MieSIrr(x,m,th);

R1 = squeeze(log10(i1));
R2 = squeeze(log10(i2));
R3 = squeeze(log10(i1+i2));

figure(1);
surf(m, x, R1);

zlabel('log_{10}(I_s / I_i)');
ylabel('Size parameter x');
xlabel('Ratio of refractive indices m');
title(strcat('Perpendicular scattered irradiance, \theta=', num2str(theta,3), ' deg'));
view([22, 28]) % useful viewing angle
shading interp  % gouraud shading of colours
lighting phong  % nice
% axis equal      % equal scales on axes
axis on         % draw axes
%light('Position',[-1 -2 0.6]);
light('Position',[0.8 -5 10]);

figure(2);
surf(m, x, R2);

zlabel('log_{10}(I_s / I_i)');
ylabel('Size parameter x');
xlabel('Ratio of refractive indices m');
title(strcat('Parallel scattered irradiance, \theta=', num2str(theta,3), ' deg'));
view([22, 28]) % useful viewing angle
shading interp  % gouraud shading of colours
lighting phong  % nice
% axis equal      % equal scales on axes
axis on         % draw axes
%light('Position',[130 -6 3]);
light('Position',[0.8 -5 10]);


%return;

figure(3);
surf(m, x, R3);
zlabel('log_{10}(I_s / I_i)');
ylabel('Size parameter x');
xlabel('Ratio of refractive indices m');
title(strcat('Total scattered irradiance, \theta=', num2str(theta,3), ' deg'));
view([22, 28]) % useful viewing angle
shading interp  % gouraud shading of colours
lighting phong  % nice
% axis equal      % equal scales on axes
axis on         % draw axes
%light('Position',[-1 -2 0.6]);
light('Position',[0.8 -5 10]);

%return;

figure(4);
surf(m, x, P);
ylabel('Size parameter x');
xlabel('Ratio of refractive indices m');
title(strcat('Polarization, \theta=', num2str(theta,3), ' deg'));
view([22, 28]) % useful viewing angle
shading interp  % gouraud shading of colours
lighting phong  % nice
% axis equal      % equal scales on axes
axis on         % draw axes
%light('Position',[-1 -2 0.6]);
light('Position',[1.5 10 10]);
