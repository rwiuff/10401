close all;
clear all;

mkdir('../MatlabFigures', 'Interferometer');

r_p = 0.055; % Reactor port opening [m]
w_0 = 0.0275; % Initial Beam Waist [m]
r = 0.125; % Minor tokamak radius [m]
freq = 60; % Frequency of probe beam [GHz]

d_0 = 0.20; % Distance between source and first lens
d_r = 0.10; % Distance between first lens and reactor wall
f_0 = 0.25; % Focal length of first lens
f_1 = 0.20; % Focal length of second lens

topstart = 2 * asin(r_p/(2 * r)) / 2;
topend = (2 * pi - 2 * topstart) / 2;
bottomstart = topend + topstart * 2;
bottomend = 2 * pi - topstart;
c = 299792458;
lambda = c / (freq * 10^9);
d = f_0 + f_1;
w_1 = (lambda * f_0) / (pi * w_0);
d_1 = (((d_0 / f_0) - 1) / ((w_0^(2) * pi) / (f_0 * lambda) + ((d_0 / f_0) - 1)^(2)) + 1) * f_0;
w_2 = (f_1 / f_0) * w_0;
d_2 = d - d_1;
d_3 = f_1 / f_0 * (f_0 + f_1 - (f_1 / f_0) * d_0);
x_0 = linspace(0, d_0, 1000);
y_0 = [];
for i = 1:1000
    w = w_0 * sqrt(1+((lambda * x_0(i)) / (pi * w_0^2))^2);
    y_0 = horzcat(y_0, w);
end

x_1 = linspace(0, d_1, 1000);
y_1 = [];
for i = 1:1000
    w = w_1 * sqrt(1+((lambda * x_1(i)) / (pi * w_1^2))^2);
    y_1 = horzcat(y_1, w);
end
y_1 = fliplr(y_1);
y = horzcat(y_0, y_1);

x_2 = linspace(0, d_2, 1000);
y_2 = [];
for i = 1:1000
    w = w_1 * sqrt(1+((lambda * x_2(i)) / (pi * w_1^2))^2);
    y_2 = horzcat(y_2, w);
end

y = horzcat(y, y_2);

x_3 = linspace(0, d_3, 1000);
y_3 = [];
for i = 1:1000
    w = w_2 * sqrt(1+((lambda * x_3(i)) / (pi * w_2^2))^2);
    y_3 = horzcat(y_3, w);
end
y_3 = fliplr(y_3);
y = horzcat(y, y_3);
for i = 1:1000
    x_1(i) = x_1(i) + d_0;
end
for i = 1:1000
    x_2(i) = x_2(i) + d_0 + d_1;
end
for i = 1:1000
    x_3(i) = x_3(i) + d_0 + d_1 + d_2;
end
x = horzcat(x_0, x_1, x_2, x_3);
testx = horzcat(x_0, x_1, x_2, x_3);
Opy = -y;
l = figure;
hold on
axis equal
p1 = plot(x, y);
p2 = plot(x, Opy, 'Color', 'red');
%p3 = xline(d_0, '--');
%xline(d_0+d_1+d_2, '--');
axPos = get(gca, 'Position');
xMinMax = xlim;
yMinMax = ylim;
zAnn = axPos(1) + ((0 - xMinMax(1)) / (xMinMax(2) - xMinMax(1))) * axPos(3);
d0Ann = axPos(1) + ((d_0 - xMinMax(1)) / (xMinMax(2) - xMinMax(1))) * axPos(3);
d1Ann = axPos(1) + ((d_0 + d_1 - xMinMax(1)) / (xMinMax(2) - xMinMax(1))) * axPos(3);
d2Ann = axPos(1) + ((d_0 + d_1 + d_2 - xMinMax(1)) / (xMinMax(2) - xMinMax(1))) * axPos(3);
d3Ann = axPos(1) + ((d_0 + d_1 + d_2 + d_3 - xMinMax(1)) / (xMinMax(2) - xMinMax(1))) * axPos(3);
yAnn = axPos(2) + ((0 - yMinMax(1)) / (yMinMax(2) - yMinMax(1))) * axPos(4);
annotation('doublearrow', [zAnn, d0Ann], [yAnn - 0.2, yAnn - 0.2]);
annotation('doublearrow', [d0Ann, d1Ann], [yAnn - 0.2, yAnn - 0.2]);
annotation('doublearrow', [d1Ann, d2Ann], [yAnn - 0.2, yAnn - 0.2]);
annotation('doublearrow', [d2Ann, d3Ann], [yAnn - 0.2, yAnn - 0.2]);
annotation('ellipse', [d0Ann - 0.01, yAnn - y_0(1000) * 1.5, 0.02, y_0(1000) * 3], 'FaceColor', 'cyan');
annotation('ellipse', [d2Ann - 0.01, yAnn - y_0(1000) * 1.5, 0.02, y_0(1000) * 3], 'FaceColor', 'cyan');
text(d_0/4, -0.22, sprintf('d_0 = %0.2f', d_0))
text(d_0+(d_1 / 4), -0.22, sprintf('d_1 = %0.2f', d_1))
text(d_0+d_1+(d_2 / 4), -0.22, sprintf('d_2 = %0.2f', d_2))
text(d_0+d_1+d_2+(d_3 / 4), -0.22, sprintf('d_3 = %0.2f', d_3))
text(d_0-0.06, 0.06, '1^{st} Lens')
text(d_0-0.06, -0.06, sprintf('f_0 = %0.2f', f_0))
text(d_0+d_1+d_2-0.06, 0.06, '2^{nd} Lens')
text(d_0+d_1+d_2-0.06, -0.06, sprintf('f_1 = %0.2f', f_1))
w1x = (d_0 + d_1);
w1x = repelem(w1x, 100);
w1y = linspace(0, w_1);
p4 = plot(w1x, w1y);
%plot(w1x(1),w1y(100), 'Marker', '^', 'MarkerSize', 2, 'MarkerEdgeColor', 'green', 'MarkerFaceColor', 'green');
%plot(w1x(1),w1y(1), 'Marker', 'V', 'MarkerSize', 2, 'MarkerEdgeColor', 'green', 'MarkerFaceColor', 'green');
yline(0, '-');
xlabel('Distance from source [m]');
ylabel('Beam Waist [m]');
title('Beam Waist In A Gaussian Beam Interferometer');
axis equal
top = linspace(topstart, topend, 100); % Creates 1000 datapoints between 0 and 4 pi
x_t = d_0 + d_r + r + r * cos(top); % x-parameter of circle
y_t = r * sin(top); % y-parameter of circle
bottom = linspace(bottomstart, bottomend, 100);
x_b = d_0 + d_r + r + r * cos(bottom); % x-parameter of circle
y_b = r * sin(bottom); % y-parameter of circle
p5 = plot(x_t, y_t, 'color', [0.911, 0.4100, 0.1700]); % Plot circle
plot(x_b, y_b, 'color', [0.911, 0.4100, 0.1700]);
lgd = legend([p1, p2, p4, p5], {'Beam Waist +', 'Beam Waist -', 'w_1', 'Reactor'});
legend('boxoff')
legend('Location', 'northwest')
lgd.NumColumns = 2;
hold off
mkdir('../MatlabFigures', 'Interferometer');
epsfilename = 'Interferometer.eps';
foldername = sprintf('../MatlabFigures/Interferometer');
fullfilename = fullfile(foldername, epsfilename);
saveas(l, fullfilename, 'epsc')
hold off
