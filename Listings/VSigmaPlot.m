clear all
close all

foldername = '../Data';
Hydrogenfile = 'Hydrogen.txt';
fullfilename = fullfile(foldername, Hydrogenfile);
H = dlmread(fullfilename);

Deuteriumfile = 'Deuterium.txt';
fullfilename = fullfile(foldername, Deuteriumfile);
D = dlmread(fullfilename);

fittype = 'power1';
str = sprintf('The type of fit is set to "%s"', fittype);
disp(str)

q = figure;
fH = fit(H(:,1), H(:,2), fittype);
pH = scatter(H(:,1), H(:,2), 'r');
hold on
pfH = plot(fH, 'r');
fD = fit(D(:,1), D(:,2), fittype);
pD = scatter(D(:,1), D(:,2), 'filled', 'b');
pfD = plot(fD, 'b');
hold off
xlabel('Voltage [kV]');
ylabel('Width of spectral line [nm]');
ytickformat('%.2f');
lgd = legend([pH, pfH, pD, pfD], ...
    {"Hydrogen Data", "Hydrogen Fit", "Deuterium Data",...
    "Deuterium Fit"});% Figure legend
legend('boxoff')
legend('Location', 'northwest')

CH = char(coeffnames(fH));
VH = coeffvalues(fH);
FH = formula(fH);

CD = char(coeffnames(fH));
VD = coeffvalues(fH);
FD = formula(fH);

for i=1:length(CH)
    FH = replace(FH, CH(i), sprintf('{%s}', num2str(VH(i))));
end

for i=1:length(CD)
    FD = replace(FD, CD(i), sprintf('{%s}', num2str(VD(i))));
end

FH = replace(FH, '*', '\cdot ');
FD = replace(FD, '*', '\cdot ');

Dollar = '$$y = ';
Dollarend = '$$';

FH = strcat(Dollar, FH, Dollarend);
FD = strcat(Dollar, FD, Dollarend);

textposH = fH(30);
textposD = fD(30);

text(30,textposH,FH, 'Interpreter', 'latex', 'FontSmoothin',...
    'on', 'HorizontalAlignment','right', 'VerticalAlignment', 'baseline');
text(30,textposD,FD, 'Interpreter', 'latex', 'FontSmoothin',...
    'on', 'HorizontalAlignment','right', 'VerticalAlignment', 'baseline');

grid on
grid minor

mkdir('../MatlabFigures', 'Asign3')
foldername = '../MatlabFigures/Asign3';
epsfilename = 'VSigma';
fullfilename = fullfile(foldername, epsfilename);
saveas(q, fullfilename, 'epsc')
str = 'Plot saved';
disp(str);
