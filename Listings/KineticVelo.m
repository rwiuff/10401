close all
clear all

foldername = '../Data';
Hydrogenfile = 'Hydrogen.txt';
fullfilename = fullfile(foldername, Hydrogenfile);
H = dlmread(fullfilename);

Deuteriumfile = 'Deuterium.txt';
fullfilename = fullfile(foldername, Deuteriumfile);
D = dlmread(fullfilename);

c = 299792458; %[m/s]

vH = [];
vD = [];
for i=1:length(H)
    vH = vertcat(vH, ((H(i,2)*10^9) * c) / (656*10^9));
end
for i=1:length(D)
    vD = vertcat(vD, ((D(i,2)*10^9) * c) / (656*10^9));
end

EKH = [];
EKD = [];
for i=1:length(H)
    EKH = vertcat(EKH, 6.242*10^(18)*(1/2)*1.6738*10^(-27)*(vH(i)^2));
end
for i=1:length(D)
    EKD = vertcat(EKD, 6.242*10^(18)*(1/2)*(2*1.6738*10^(-27))*(vD(i)^2));
end

q1 = figure;
yyaxis left
eqn = fittype('a*sqrt(x)+b');
str = sprintf('The type of fit is set to "%s"', eqn);
disp(str)
hold on
vHFit = fit(H(:,1), vH, eqn, 'StartPoint', [0,0]);
vDFit = fit(D(:,1), vD, eqn, 'StartPoint', [0,0]);
pSH = scatter(H(:,1), vH, 'filled');
pSD = scatter(D(:,1), vD);
pVH = plot(vHFit, 'b');
pVD = plot(vDFit, '-.');
ylabel('Velocity [m/s]');
yyaxis right
eqn = fittype('a*x+b');
str = sprintf('The type of fit is set to "%s"', eqn);
disp(str)
EKHFit = fit(H(:,1), EKH, eqn, 'StartPoint', [0,0]);
EKDFit = fit(D(:,1), EKD, eqn, 'StartPoint', [0,0]);
pEH = scatter(H(:,1), EKH, 'filled');
pED = scatter(D(:,1), EKD);
peH = plot(EKHFit);
peD = plot(EKDFit, '-.');
hold off
xlabel('Voltage [-kV]');
ylabel('Energy [eV]');
% ylim([EKH(8)*1.01,EKH(1)*1.01])
lgd = legend([pSH, pSD, pVH, pVD, pEH, pED, peH, peD], ...
     {"Hydrogen velocities", "Deuterium velocities", "Hydrogen velocity fit",...
     "Deuterium velocity fit", "Hydrogen energies", "Deuterium energies",...
     "Hydrogen energy fit", "Deuterium energy fit"});% Figure legend
legend('boxoff')
legend('Location', 'southoutside')
lgd.NumColumns = 2;
grid on
grid minor

mkdir('../MatlabFigures', 'Asign3')
foldername = '../MatlabFigures/Asign3';
epsfilename = 'KineticVelo';
fullfilename = fullfile(foldername, epsfilename);
saveas(q1, fullfilename, 'epsc')
str = 'Plot saved';
disp(str);
disp(vHFit)
disp(vDFit)
disp(EKHFit)
disp(EKDFit)
