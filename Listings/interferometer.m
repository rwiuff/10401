close;
clear;

w_0 = 0.0275; % Initial Beam Waist [m]
r = 0.125; % Minor tokamak radius [m]
freq = 60; % Frequency of probe beam [GHz]

d_1 = 0.01; % Distance between source and first lens
d_2 = 0.01; % Distance between first and second lens
d_3 = 0.01; % Distance between second lens and chamber end
f_1 = 0.01; % Focal length of first lens
f_2 = 0.01; % Focal length of second lens

c = 299792458; % The speed of light [m/s]
lambda = c/(freq*10^9); % Wavelength [m]

x = linspace(0,1,1000);
y = [];
for i=0:1000
    y(i) = w_0*sqrt(1+((lambda*x(i))/(pi*w_0^2))^2);
end
plot(x,y)


t=linspace(0,4*pi, 1000); % Creates 1000 datapoints between 0 and 4 pi
x_1 = r*cos(t); % x-parameter of circle
y_1 = r*sin(t); % y-parameter of circle
q = figure; % Figure variable
hold on
plot(x_1,y_1); % Plot circle
xlabel('Distance [m]');
ylabel('Distance [m]');
title('Profile of Gaussian Beam Interferometer');
legend('Tokamak chamber wall');
axis equal
hold off
mkdir('../MatlabFigures', 'Interferometer');
epsfilename = 'Interferometer.eps';
foldername = sprintf('../MatlabFigures/Interferometer');
fullfilename = fullfile(foldername,epsfilename);
saveas(q, fullfilename, 'epsc')
