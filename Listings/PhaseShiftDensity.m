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

y = linspace(n_0,n_1,1000);
Phi0 = [];
Phi1 = [];
Phi2 = [];
for i=1:1000
    Phi0 = horzcat(Phi0, 8.416e-8*(y(i)/(f_0*10^9)));
end
for i=1:1000
    Phi1 = horzcat(Phi1, 8.416e-8*(y(i)/(f_1*10^9)));
end
for i=1:1000
    Phi2 = horzcat(Phi2, 8.416e-8*(y(i)/(f_2*10^9)));
end
l = figure; % Creates figure
hold on
p0 = plot(Phi0,y);
p1 = plot(Phi1,y);
p2 = plot(Phi2,y);
xlabel('Phase shift [2\pi]'); % x-axis label
ylabel('Electron Density [m^-3]'); % y-axis label
%title('Electron Density as function of measured phase shift'); % Figure title
lgd = legend([p0, p1, p2], ...
    {'60GHz', '98GHz', '130GHz'});% Figure legend
legend('boxoff')
legend('Location', 'northwest')
hold off
%-------------------------------------------------------------------------%


%------------------Saving figure as Encapsulated Postscript---------------%
epsfilename = 'PhaseShift.eps'; % Savename for the figure
foldername = sprintf('../MatlabFigures/PhaseShift'); % Folder path
fullfilename = fullfile(foldername, epsfilename); % Filename path
saveas(l, fullfilename, 'epsc') % Save the figure as eps
%-------------------------------------------------------------------------%

