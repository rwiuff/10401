close all;
clear all;

mkdir('../MatlabFigures', 'Interferometer'); % Create save directory

%--------------------------Input Parameters-------------------------------%
r_p = 0.055; % Reactor port opening [m]
w_0 = 0.0275; % Initial Beam Waist [m]
r = 0.125; % Minor tokamak radius [m]
freq = 60; % Frequency of probe beam [GHz]

d_0 = 0.20; % Distance between source and first lens
d_r = 0.10; % Distance between first lens and reactor wall
f_0 = 0.25; % Focal length of first lens
f_1 = 0.20; % Focal length of second lens
%-------------------------------------------------------------------------%

%----------Distances and Lowest Beam Waists Calculations------------------%
c = 299792458; % The speec of light [m/s]
lambda = c / (freq * 10^9); % Calculates the wavelength [m]
d = f_0 + f_1; % Calculates the distance between lenses
w_1 = (lambda * f_0) / ...
    (pi * w_0);% Calculates the beam waist between lenses
d_1 = (((d_0 / f_0) - 1) / ((w_0^(2) * pi) / ...
    (f_0 * lambda) + ((d_0 / f_0) - 1)^(2)) + 1) * ...
    f_0;% Distance after 1st lens to lowest beam waist
w_2 = (f_1 / f_0) * w_0; % Calculates w_2 after 2nd second lens
d_2 = d - d_1; % Distance from w_2 to second lens
d_3 = f_1 / f_0 *...
    (f_0 + f_1 -...
    (f_1 / f_0) * d_0);% Distance to lowest beam waist after 2nd lens
%-------------------------------------------------------------------------%

%------------------------Beam Waist Datapoints----------------------------%
x_0 = linspace(0, d_0, 1000);% Plotpoints till first lens
y_0 = [];
for i = 1:1000
    w = w_0 * ...
        sqrt(1+((lambda * x_0(i)) / ...
        (pi * w_0^2))^2);% Calculates the beam waist from source
    y_0 = horzcat(y_0, w);
end

x_1 = linspace(0, d_1, 1000); % Plotpoints from first lens to w_1
y_1 = [];
for i = 1:1000
    w = w_1 * ...
        sqrt(1+((lambda * x_1(i)) / ...
        (pi * w_1^2))^2);% Beam waist between w_1 and the 1st lens
    y_1 = horzcat(y_1, w);
end
y_1 = fliplr(y_1); % Flips the y-plot points as the beam waist is declining
y = horzcat(y_0, y_1);

x_2 = linspace(0, d_2, 1000); % Plotpoints from w_1 to second lens
y_2 = [];
for i = 1:1000
    w = w_1 * ...
        sqrt(1+((lambda * x_2(i)) / ...
        (pi * w_1^2))^2);% Beam waist between w_1 and the 2nd lens
    y_2 = horzcat(y_2, w);
end

y = horzcat(y, y_2);

x_3 = linspace(0, d_3, 1000); % Plotpoints from 2nd lens to w_2
y_3 = [];
for i = 1:1000
    w = w_2 * ...
        sqrt(1+((lambda * x_3(i)) / ...
        (pi * w_2^2))^2);% Evolving beam waist after the 2nd lens to w_2
    y_3 = horzcat(y_3, w);
end
y_3 = fliplr(y_3); % Flips the y-plot points as the beam waist is declining
y = horzcat(y, y_3);

for i = 1:1000
    x_1(i) = x_1(i) + d_0; % Creates x-axis data points
end
for i = 1:1000
    x_2(i) = x_2(i) + d_0 + d_1; % Creates x-axis data points
end
for i = 1:1000
    x_3(i) = x_3(i) + d_0 + d_1 + d_2; % Creates x-axis data points
end
x = horzcat(x_0, x_1, x_2, x_3);

Opy = -y; % Flip beam waist to plot 'Beam Waist -'
%-------------------------------------------------------------------------%

%---------------------Beam Waist Propagation Plot-------------------------%
l = figure; % Creates figure
hold on
axis equal
p1 = plot(x, y); % Plots Beam Waist+
p2 = plot(x, Opy, 'Color', 'red'); % Plots Beam Waist-
%p3 = xline(d_0, '--');
%xline(d_0+d_1+d_2, '--');
axPos = get(gca, 'Position'); % Get normalised axis coordinates
xMinMax = xlim; % Get normalised axis coordinates
yMinMax = ylim; % Get normalised axis coordinates
zAnn = axPos(1) + ((0 - xMinMax(1)) / (xMinMax(2) - xMinMax(1))) ...
    * axPos(3);% Defines points relative to axis in normalised coordinates
d0Ann = axPos(1) + ((d_0 - xMinMax(1)) / (xMinMax(2) - xMinMax(1))) ...
    * axPos(3);% Defines points relative to axis in normalised coordinates
d1Ann = axPos(1) + ((d_0 + d_1 - xMinMax(1)) / (xMinMax(2) - xMinMax(1))) ...
    * axPos(3);% Defines points relative to axis in normalised coordinates
d2Ann = axPos(1) + ((d_0 + d_1 + d_2 - xMinMax(1)) / (xMinMax(2) - xMinMax(1))) ...
    * axPos(3);% Defines points relative to axis in normalised coordinates
d3Ann = axPos(1) + ((d_0 + d_1 + d_2 + d_3 - xMinMax(1)) / (xMinMax(2) - xMinMax(1))) ...
    * axPos(3);% Defines points relative to axis in normalised coordinates
yAnn = axPos(2) + ((0 - yMinMax(1)) / (yMinMax(2) - yMinMax(1))) ...
    * axPos(4);% Defines points relative to axis in normalised coordinates
annotation('doublearrow', ...
    [zAnn, d0Ann], [yAnn - 0.2, yAnn - 0.2]);% Annotates d_0
annotation('doublearrow', ...
    [d0Ann, d1Ann], [yAnn - 0.2, yAnn - 0.2]);% Annotates d_1
annotation('doublearrow', ...
    [d1Ann, d2Ann], [yAnn - 0.2, yAnn - 0.2]);% Annotates d_2
annotation('doublearrow', ...
    [d2Ann, d3Ann], [yAnn - 0.2, yAnn - 0.2]);% Annotates d_3
annotation('ellipse', ...
    [d0Ann - 0.01, yAnn - y_0(1000) * 1.5, 0.02, y_0(1000) * 3], ...
    'FaceColor', 'cyan');% Draws 1st lens
annotation('ellipse', ...
    [d2Ann - 0.01, yAnn - y_0(1000) * 1.5, 0.02, y_0(1000) * 3], ...
    'FaceColor', 'cyan');% Draws 2nd lens
text(d_0/4, -0.22, sprintf('d_0 = %0.2f', d_0)) % Distance annotation
text(d_0+(d_1 / 4), -0.22, ...
    sprintf('d_1 = %0.2f', d_1))% Distance annotation
text(d_0+d_1+(d_2 / 4), -0.22, ...
    sprintf('d_2 = %0.2f', d_2))% Distance annotation
text(d_0+d_1+d_2+(d_3 / 4), -0.22, ...
    sprintf('d_3 = %0.2f', d_3))% Distance annotation
text(d_0-0.06, 0.06, '1^{st} Lens') % Lens annotation
text(d_0-0.06, -0.06, sprintf('f_0 = %0.2f', f_0)) % Lens annotation
text(d_0+d_1+d_2-0.06, 0.06, '2^{nd} Lens') % Lens annotation
text(d_0+d_1+d_2-0.06, -0.06, ...
    sprintf('f_1 = %0.2f', f_1))% Lens annotation
w1x = (d_0 + d_1); % w_1 coordinate
w1x = repelem(w1x, 100); % w_1 x-coordinate array
w1y = linspace(0, w_1); % w_1 y-coordinate array
p4 = plot(w1x, w1y); % Plots w_1
yline(0, '-'); % Plots optical axis
xlabel('Distance from source [m]'); % x-axis label
ylabel('Beam Waist [m]'); % y-axis label
title('Beam Waist In A Gaussian Beam Interferometer'); % Figure title

topstart = 2 * asin(r_p/(2 * r)) / 2; % Calculates reactor plot coordinates
topend = (2 * pi - 2 * topstart) / 2; % Calculates reactor plot coordinates
bottomstart = topend + topstart * 2; % Calculates reactor plot coordinates
bottomend = 2 * pi - topstart; % Calculates reactor plot coordinates

top = linspace(topstart, ...
    topend, 100);% Creates 100 datapoints for the top half reactor cutout
x_t = d_0 + d_r + r + r * cos(top); % x-parameter of circle
y_t = r * sin(top); % y-parameter of circle
bottom = linspace(bottomstart, ...
    bottomend, 100);% Creates 100 datapoints for the bottom half reactor
x_b = d_0 + d_r + r + r * cos(bottom); % x-parameter of circle
y_b = r * sin(bottom); % y-parameter of circle
p5 = plot(x_t, y_t, 'color', ...
    [0.911, 0.4100, 0.1700]);% Plot top half reactor
plot(x_b, y_b, 'color', ...
    [0.911, 0.4100, 0.1700]);% Plot bottom half reactor
lgd = legend([p1, p2, p4, p5], ...
    {'Beam Waist +', 'Beam Waist -', 'w_1', 'Reactor'});% Figure legend
legend('boxoff')
legend('Location', 'northwest')
lgd.NumColumns = 2;
hold off
%-------------------------------------------------------------------------%

%------------------Saving figure as Encapsulated Postscript---------------%
epsfilename = 'Interferometer.eps'; % Savename for the figure
foldername = sprintf('../MatlabFigures/Interferometer'); % Folder path
fullfilename = fullfile(foldername, epsfilename); % Filename path
saveas(l, fullfilename, 'epsc') % Save the figure as eps
%-------------------------------------------------------------------------%

%----------------------------Print Outputs--------------------------------%
w_1
w_2
d_1
d_2
d_3
%-------------------------------------------------------------------------%
