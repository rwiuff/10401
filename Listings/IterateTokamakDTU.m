close;
clear;

titl = ["Desired output power [MW]", "Maximum wall load [MW m^-2]",...
    "Magnetic field at the edge of the coil [T]",...
    "Tensile strenght of the magnetic field coils [atm]"];
l = 5;
p1 = [];
x = [];
mkdir('../Matlab Figs', titl(l))

for i = 2000:5000
    [b, c, a, R_0, A, A_p, V_p, P_dens, p,...
        n, B_0, beta, tau_E_min, C_per_watt] =...
        tokamakDTU_asign_1(0.01, 1, 2, 1000, 4, 13, 3000, 0.4);
    q=[b, c, a, R_0, A, A_p, V_p, P_dens, p,...
        n, B_0, beta, tau_E_min, C_per_watt];
    p1=cat(1,p1,q);
    x=cat(1,x,i);
end

T = ["Blanket-shield thickness [m]", "Magnet coil thickness [m]"...
    "Minor radius [m]", "Major radius [m]", "Aspect ratio []"...
    "Plasma surface [m^2]", "Plasma volume [m^3]", "Power density [W m^-1]"...
    "Plasma pressure [Pa]", "Particle density [m^-3]",...
    "Magnetic field at magnetic axis [T]", "Plasma beta in the centre []"...
    "Min confinement time for satisfaction of (p tau_E)_min [s]",...
    "The cost of the powerplant [$]"];
for k = 1:14
    q = figure;
    y = p1(:,k);
    CM = jet(14);
    plot(x, y,'color', CM(k,:));
    ylabel(T(k));
    xlabel(titl(l));
    ytickformat('%.2f');
    epsfilename = sprintf('%s.eps',T(k));
    foldername = sprintf('../Matlab Figs/%s', titl(l));
    fullfilename = fullfile(foldername,epsfilename);
    saveas(q, fullfilename, 'epsc')
end

