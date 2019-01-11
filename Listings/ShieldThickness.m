close;
clear;

mkdir('../Matlab Figs', 'ShieldThickness');

kappa = 2; % Ratio between semi axes
a_min = 1; % Smallest semi axis
b = 1.2; % Wall thickness at the axis
t=linspace(0,4*pi, 1000); % Creates 1000 datapoints between 0 and 4 pi
x_1 = a_min*cos(t); % x-parameter of small ellipse
y_1 = kappa*a_min*sin(t); % y-parameter of small ellipse
x_2 = (b+a_min)*cos(t); % x-parameter of large ellipse
y_2 = kappa*(b+a_min)*sin(t); % y-parameter of large ellipse
q = figure; % Designate figure
subplot(1,2,1); % Create subplot
hold on
plot(x_1,y_1); % Plot small ellipse
plot(x_2,y_2); % Plot large ellipse
xlabel('Profile of the blanket/shield');
legend('Inner wall','Outer wall');
axis equal
x = []; % x-component array
y = []; % y-component array
f = []; % Wall thickness array
for i=1:1000 % Iteration over thicknesses
    x(i) = x_2(i) - x_1(i); % Finds the difference in x-components
    y(i) = y_2(i) - y_1(i); % Finds the difference in y-components
    f(i) = sqrt(x(i)^2+y(i)^2); % Finds the distances between x,y points
end
subplot(1,2,2); % Create subplot
plot(t,f) % Plot distances
xlabel('Degrees around the plasma')
ylabel('Thickness of the shield')
grid on, axis ([0 2*pi 0 3])
hold off
epsfilename = 'ShieldThickness.eps';
foldername = sprintf('../Matlab Figs/ShieldThickness');
fullfilename = fullfile(foldername,epsfilename);
saveas(q, fullfilename, 'epsc')
