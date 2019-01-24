close all
clear all

foldername = '../Data';
Hydrogenfile = 'Contours.txt';
fullfilename = fullfile(foldername, Hydrogenfile);
C = dlmread(fullfilename);

I = C(:,1);
V = C(:,2);
P = C(:,3);
Cts = C(:,4);
Flux = C(:,5);
ns = C(:,6);

x = I;
y = V;
z = Cts;
[X,Y]=meshgrid(x,y);
Z = reshape(z,numel(y),numel(x));
contourf(X,Y,Z)

% mkdir('../MatlabFigures', 'Asign3')
% foldername = '../MatlabFigures/Asign3';
% epsfilename = 'PlasmaUI.eps';
% fullfilename = fullfile(foldername, epsfilename);
% saveas(q1, fullfilename, 'epsc')
% str = 'Plot saved';
% disp(str);
