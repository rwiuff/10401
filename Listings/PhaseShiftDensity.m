close all;
clear all;

mkdir('../MatlabFigures', 'PhaseShift'); % Create save directory

%--------------------------Input Parameters-------------------------------%
f_0 = 60; % 1st beam frequency [GHz]
f_1 = 98; % 2nd beam frequency [GHz]
f_2 = 130; % 3rd beam frequency [GHz]

n_0 = 10^16; % Lower electron density
n_1 = 10^18; % Hihger electron density
%-------------------------------------------------------------------------%

c = 299792458; % The speec of light [m/s]

y = linspace(n_0, n_1, 1000); % Datapoints between n_0 and n_1
Phi0 = []; % Phase shift array
Phi1 = []; % Phase shift array
Phi2 = []; % Phase shift array

for i = 1:1000 % Calculates Phase Shifts for frequency = f_0
    Phi0 = horzcat(Phi0, 8.416e-8*(y(i) / (f_0 * 10^9)));
end
for i = 1:1000 % Calculates Phase Shifts for frequency = f_1
    Phi1 = horzcat(Phi1, 8.416e-8*(y(i) / (f_1 * 10^9)));
end
for i = 1:1000 % Calculates Phase Shifts for frequency = f_2
    Phi2 = horzcat(Phi2, 8.416e-8*(y(i) / (f_2 * 10^9)));
end

l = figure; % Creates figure
hold on
p0 = plot(Phi0, y); % Plot electron denisty as function of phase shift
p1 = plot(Phi1, y); % Plot electron denisty as function of phase shift
p2 = plot(Phi2, y); % Plot electron denisty as function of phase shift
xlabel('Phase shift [^{1}/_{2\pi}]'); % x-axis label
ylabel('Electron Density [m^-3]'); % y-axis label
%title('Electron Density as function of measured phase shift'); % Figure title
lgd = legend([p0, p1, p2], ...
    {'60GHz', '98GHz', '130GHz'});% Figure legend
legend('boxoff')
legend('Location', 'northwest')
grid on
grid minor
hold off

%------------------Saving figure as Encapsulated Postscript---------------%
epsfilename = 'PhaseShift.eps'; % Savename for the figure
foldername = sprintf('../MatlabFigures/PhaseShift'); % Folder path
fullfilename = fullfile(foldername, epsfilename); % Filename path
saveas(l, fullfilename, 'epsc') % Save the figure as eps
%-------------------------------------------------------------------------%
