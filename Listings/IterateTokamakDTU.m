close all
clear all

titl = ["Desired output power [MW]", "Maximum wall load [MW m^-2]",...
    "Magnetic field at the edge of the coil [T]",...
    "Tensile strenght of the magnetic field coils [atm]"];
foldertitl = "TokamakIterations";
l = 1; % Choose titl
p1 = [];
x = [];
mkdir('../MatlabFigures', foldertitl)

if l == 1
    x = linspace(1,3000,1000);
    for i = 1:1000
    [b, c, a, R_0, A, A_p, V_p, P_dens, p,...
        n, B_0, beta, tau_E_min, C_per_watt] =...
        tokamakDTU_asign_1(0.01, 1, 2, x(i), 4, 13, 3000, 0.4);
    q=[b, c, a, R_0, A, A_p, V_p, P_dens, p,...
        n, B_0, beta, tau_E_min, C_per_watt];
    p1=cat(1,p1,q);
    end
elseif l == 2
    for i = 1:1000:10
    [b, c, a, R_0, A, A_p, V_p, P_dens, p,...
        n, B_0, beta, tau_E_min, C_per_watt] =...
        tokamakDTU_asign_1(0.01, 1, 2, 1000, i, 13, 3000, 0.4);
    q=[b, c, a, R_0, A, A_p, V_p, P_dens, p,...
        n, B_0, beta, tau_E_min, C_per_watt];
    p1=cat(1,p1,q);
    x=cat(1,x,i);
    end
elseif l == 3
    for i = 10:1000:20
    [b, c, a, R_0, A, A_p, V_p, P_dens, p,...
        n, B_0, beta, tau_E_min, C_per_watt] =...
        tokamakDTU_asign_1(0.01, 1, 2, 1000, 4, i, 3000, 0.4);
    q=[b, c, a, R_0, A, A_p, V_p, P_dens, p,...
        n, B_0, beta, tau_E_min, C_per_watt];
    p1=cat(1,p1,q);
    x=cat(1,x,i);
    end
elseif l == 4
    for i = 2000:1000:5000
    [b, c, a, R_0, A, A_p, V_p, P_dens, p,...
        n, B_0, beta, tau_E_min, C_per_watt] =...
        tokamakDTU_asign_1(0.01, 1, 2, 1000, 4, 13, i, 0.4);
    q=[b, c, a, R_0, A, A_p, V_p, P_dens, p,...
        n, B_0, beta, tau_E_min, C_per_watt];
    p1=cat(1,p1,q);
    x=cat(1,x,i);
    end
end
b = p1(:,1);
c = p1(:,2);
a = p1(:,3);
R_0 = p1(:,4);
A = p1(:,5);
A_p = p1(:,6);
V_p = p1(:,7);
P_dens = p1(:,8);
p = p1(:,9);
n = p1(:,10);
B_0 = p1(:,11);
beta = p1(:,12);
tau_E_min = p1(:,13);
C_per_watt = p1(:,14);


fileT = ["Blanket-shieldThickness", "MagnetCoilThickness"...
    "MinorRadius", "MajorRadius", "AspectRatio"...
    "PlasmaSurface", "PlasmaVolume", "PowerDensity"...
    "PlasmaPressure", "ParticleDensity",...
    "MagneticFieldAtMagneticAxis", "PlasmaBetaInTheCentre"...
    "MinConfinementTime",...
    "TheCostOfThePowerplant"];

for k = 1:14
    q = figure;
    y = p1(:,k);
    CM = jet(14);
    plot(x, y,'color', CM(k,:));
    xlabel(titl(l));
    ytickformat('%.2f');
    epsfilename = sprintf('%s.eps',fileT(k));
    foldername = sprintf('../MatlabFigures/%s', foldertitl);
    fullfilename = fullfile(foldername,epsfilename);
    saveas(q, fullfilename, 'epsc')
end

lgd = legend([p1, p2, p4, p5], ...
    {"Blanket-shield thickness [m]", "Magnet coil thickness [m]"...
    "Minor radius [m]", "Major radius [m]", "Aspect ratio []"...
    "Plasma surface [m^2]", "Plasma volume [m^3]", "Power density [W m^-1]"...
    "Plasma pressure [Pa]", "Particle density [m^-3]",...
    "Magnetic field at magnetic axis [T]", "Plasma beta in the centre []"...
    "Min confinement time for satisfaction of (p tau_E)_min [s]",...
    "The cost of the powerplant [$]"});% Figure legend
legend('boxoff')
legend('Location', 'northwest')

