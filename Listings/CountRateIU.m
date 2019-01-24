close all
clear all

foldername = '../Data';
Hydrogenfile = 'Hydrogen.txt';
fullfilename = fullfile(foldername, Hydrogenfile);
H = dlmread(fullfilename);

Deuteriumfile = 'Deuterium.txt';
fullfilename = fullfile(foldername, Deuteriumfile);
D = dlmread(fullfilename);




mkdir('../MatlabFigures', 'Asign3')
foldername = '../MatlabFigures/Asign3';
epsfilename = 'PlasmaUI.eps';
fullfilename = fullfile(foldername, epsfilename);
saveas(q1, fullfilename, 'epsc')
str = 'Plot saved';
disp(str);
