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
Q = C(:,7);

x = I;
y = V;
z = Cts;
Tickcount = [];
Ticklabel = [];
Tarray = sort(z);
for i=1:length(z)
    Tickcount = [Tickcount, Tarray(i)];
    Ticklabel = [Ticklabel, sprintf("%0.2f", Tarray(i))];
end
q = figure;
hold on
F = scatteredInterpolant(x, y, z, 'nearest');
[xq,yq] = meshgrid(min(x)-1:0.01:max(x)+2,min(y)-1:0.01:max(y)+1);
cq = F(xq,yq);
h = pcolor(xq,yq,cq);
h.EdgeColor = 'none';
scatter(x',y',50,z','filled',...
    'MarkerFaceAlpha', .75, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5)
c = colorbar('eastoutside', 'Ticks', Tickcount, 'TickLabels', Ticklabel);
ylabel('Applied voltage [-kV]')
xlabel('Current [mA]')
c.Label.String = 'Neutron count rates [counts/s]';
hold off

mkdir('../MatlabFigures', 'Asign3')
foldername = '../MatlabFigures/Asign3';
epsfilename = 'Counts.eps';
fullfilename = fullfile(foldername, epsfilename);
saveas(q, fullfilename, 'epsc')
str = 'Plot 1 of 4 saved';
disp(str);
close all

x = I;
y = V;
z = Flux;
Tickcount = [];
Ticklabel = [];
Tarray = sort(z);
for i=1:length(z)
    Tickcount = [Tickcount, Tarray(i)];
    Ticklabel = [Ticklabel, sprintf("%0.2e", Tarray(i))];
end
q = figure;
hold on
F = scatteredInterpolant(x, y, z, 'nearest');
[xq,yq] = meshgrid(min(x)-1:0.01:max(x)+2,min(y)-1:0.01:max(y)+1);
cq = F(xq,yq);
h = pcolor(xq,yq,cq);
h.EdgeColor = 'none';
scatter(x',y',50,z','filled',...
    'MarkerFaceAlpha', .75, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5)
c = colorbar('eastoutside', 'Ticks', Tickcount, 'TickLabels', Ticklabel);
ylabel('Applied voltage [-kV]')
xlabel('Current [mA]')
c.Label.String = 'Neutron flux at detector [s^{-1}m^{-2}]';
decrease_by = 0.01;
axpos = get(gca,'position');
axpos(3) = axpos(3) - decrease_by;
set(gca,'position',axpos);

hold off

mkdir('../MatlabFigures', 'Asign3')
foldername = '../MatlabFigures/Asign3';
epsfilename = 'NFlux.eps';
fullfilename = fullfile(foldername, epsfilename);
saveas(q, fullfilename, 'epsc')
str = 'Plot 2 of 4 saved';
disp(str);
close all

x = I;
y = V;
z = P;
Tickcount = [];
Ticklabel = [];
Tarray = sort(z);
Tickcount = [Tarray(1), Tarray(4), Tarray(5)];
Ticklabel = [sprintf("%0.2f", Tarray(1)),...
    sprintf("%0.2f", Tarray(4)),...
    sprintf("%0.2f", Tarray(5))];
q = figure;
hold on
F = scatteredInterpolant(x, y, z, 'nearest');
[xq,yq] = meshgrid(min(x)-1:0.01:max(x)+2,min(y)-1:0.01:max(y)+1);
cq = F(xq,yq);
h = pcolor(xq,yq,cq);
h.EdgeColor = 'none';
scatter(x',y',50,z','filled',...
    'MarkerFaceAlpha', .75, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5)
c = colorbar('eastoutside', 'Ticks', Tickcount, 'TickLabels', Ticklabel);
ylabel('Applied voltage [-kV]')
xlabel('Current [mA]')
c.Label.String = 'Fusion Power [W]';
hold off

mkdir('../MatlabFigures', 'Asign3')
foldername = '../MatlabFigures/Asign3';
epsfilename = 'P.eps';
fullfilename = fullfile(foldername, epsfilename);
saveas(q, fullfilename, 'epsc')
str = 'Plot 3 of 4 saved';
disp(str);
close all


x = I;
y = V;
z = Q;
Tickcount = [];
Ticklabel = [];
Tarray = sort(z);
for i=1:length(z)
    Tickcount = [Tickcount, Tarray(i)];
    Ticklabel = [Ticklabel, sprintf("%0.2e", Tarray(i))];
end
q = figure;
hold on
F = scatteredInterpolant(x, y, z, 'nearest');
[xq,yq] = meshgrid(min(x)-1:0.01:max(x)+2,min(y)-1:0.01:max(y)+1);
cq = F(xq,yq);
h = pcolor(xq,yq,cq);
h.EdgeColor = 'none';
scatter(x',y',50,z','filled',...
    'MarkerFaceAlpha', .75, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5)
c = colorbar('eastoutside', 'Ticks', Tickcount, 'TickLabels', Ticklabel);
ylabel('Applied voltage [-kV]')
xlabel('Current [mA]')
c.Label.String = 'Q-factor';

hold off

mkdir('../MatlabFigures', 'Asign3')
foldername = '../MatlabFigures/Asign3';
epsfilename = 'Q.eps';
fullfilename = fullfile(foldername, epsfilename);
saveas(q, fullfilename, 'epsc')
str = 'Plot 4 of 4 saved';
disp(str);
close all
str = 'Done';
disp(str)
