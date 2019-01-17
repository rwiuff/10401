close all;
clear all;

mkdir('../MatlabFigures', 'Interferometer');

w_0 = 0.0275; % Initial Beam Waist [m]
r = 0.125; % Minor tokamak radius [m]
freq = 60; % Frequency of probe beam [GHz]

d_0 = 0.1; % Distance between source and first lens
f_0 = 0.2; % Focal length of first lens
f_1 = 0.05; % Focal length of second lens

c = 299792458;
lambda = c/(freq*10^9)
d = f_0+f_1;
w_1 = (lambda*f_0)/(pi*w_0)
d_1 = (((d_0/f_0)-1)/((w_0^(2)*pi)/(f_0*lambda)+((d_0/f_0)-1)^(2))+1)*f_0
w_2 = (f_1/f_0)*w_0
d_2 = d - d_1
d_3 = f_1/f_0*(f_0+f_1-(f_1/f_0)*d_0)
x_0 = linspace(0,d_0,1000);
y_0 = [];
for i = 1:1000
    w = w_0*sqrt(1+((lambda*x_0(i))/(pi*w_0^2))^2);
    y_0 = horzcat(y_0, w);
end

x_1 = linspace(0,d_1,1000);
y_1 = [];
for i = 1:1000
    w = w_1*sqrt(1+((lambda*x_1(i))/(pi*w_1^2))^2);
    y_1 = horzcat(y_1, w);
end
y_1 = fliplr(y_1);
y = horzcat(y_0,y_1);

x_2 = linspace(0,d_2,1000);
y_2 = [];
for i = 1:1000
    w = w_1*sqrt(1+((lambda*x_2(i))/(pi*w_1^2))^2);
    y_2 = horzcat(y_2, w);
end

y = horzcat(y,y_2);

x_3 = linspace(0,d_3,1000);
y_3 = [];
for i = 1:1000
    w = w_2*sqrt(1+((lambda*x_3(i))/(pi*w_2^2))^2);
    y_3 = horzcat(y_3, w);
end
y_3 = fliplr(y_3);
y = horzcat(y,y_3);
for i = 1:1000
    x_1(i) = x_1(i)+d_0;
end
for i = 1:1000
    x_2(i) = x_2(i)+d_0+d_1;
end
for i = 1:1000
    x_3(i) = x_3(i)+d_0+d_1+d_2;
end
x = horzcat(x_0,x_1,x_2,x_3);
testx = horzcat(x_0,x_1,x_2,x_3);
Opy = -y;
l = figure;
hold on
axis equal
plot(x,y)
plot(x,Opy)
xline(d_0, '--', {'First lens'});
xline(d_0+d_1+d_2, '--', {'Second lens'});
axPos = get(gca, 'Position');
xMinMax = xlim;
yMinMax = ylim;
zAnn = axPos(1) + ((0 - xMinMax(1))/(xMinMax(2)-xMinMax(1))) * axPos(3);
d0Ann = axPos(1) + ((d_0 - xMinMax(1))/(xMinMax(2)-xMinMax(1))) * axPos(3);
d1Ann = axPos(1) + ((d_0 + d_1 - xMinMax(1))/(xMinMax(2)-xMinMax(1))) * axPos(3);
d2Ann = axPos(1) + ((d_0 + d_1 + d_2 - xMinMax(1))/(xMinMax(2)-xMinMax(1))) * axPos(3);
d3Ann = axPos(1) + ((d_0 + d_1 + d_2 + d_3 - xMinMax(1))/(xMinMax(2)-xMinMax(1))) * axPos(3);
yAnn = axPos(2) + ((-0.05 - yMinMax(1))/(yMinMax(2)-yMinMax(1))) * axPos(4);
annotation('doublearrow', [zAnn, d0Ann], [yAnn,yAnn]);
annotation('doublearrow', [d0Ann, d1Ann], [yAnn,yAnn]);
annotation('doublearrow', [d1Ann, d2Ann], [yAnn,yAnn]);
annotation('doublearrow', [d2Ann, d3Ann], [yAnn,yAnn]);
xlabel('Distance [m]');
ylabel('Beam Waist [m]');
title('Profile of Gaussian Beam Interferometer');
axis equal
hold off
mkdir('../MatlabFigures', 'Interferometer');
epsfilename = 'Interferometer.eps';
foldername = sprintf('../MatlabFigures/Interferometer');
fullfilename = fullfile(foldername,epsfilename);
saveas(l, fullfilename, 'epsc')
hold off

% t=linspace(0,4*pi, 1000); % Creates 1000 datapoints between 0 and 4 pi
% x_1 = r + r*cos(t); % x-parameter of circle
% y_1 = r + r*sin(t); % y-parameter of circle
% q = figure; % Figure variable
% hold on
% plot(x_1,y_1); % Plot circle
% xlabel('Distance [m]');
% ylabel('Distance [m]');
% title('Profile of Gaussian Beam Interferometer');
% legend('Tokamak chamber wall');
% axis equal
% hold off
% mkdir('../MatlabFigures', 'Interferometer');
% epsfilename = 'Interferometer.eps';
% foldername = sprintf('../MatlabFigures/Interferometer');
% fullfilename = fullfile(foldername,epsfilename);
% saveas(q, fullfilename, 'epsc')
