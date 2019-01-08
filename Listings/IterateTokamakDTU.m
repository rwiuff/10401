%**************************************************************************
% Name:         IterateTokamakDTU
%
% Version:      1.0
%
%**************************************************************************

close;
clear;

p1 = [];
x = [];

for i = 500:1500
    [b, c, a, R_0, A, A_p, V_p, P_dens, p,...
        n, B_0, beta, tau_E_min, C_per_watt] =...
        tokamakDTU_asign_1(0.01, 1, 2, i, 4, 13, 3000, 0.4);
    q=[b, c, a, R_0, A, A_p, V_p, P_dens, p,...
        n, B_0, beta, tau_E_min, C_per_watt];
    p1=cat(1,p1,q);
    x=cat(1,x,i);
end

T = ["Blanket/shield thickness [m]", "Magnet coil thickness [m]"...
    "Minor radius [m]", "Major radius [m]", "Aspect ratio []"...
    "Plasma surface [m^2]", "Plasma volume [m^3]", "Power density [W/m]"...
    "Plasma pressure [Pa]", "Particle density [m^-3]",...
    "Magnetic field at magnetic axis [T]", "Plasma beta in the centre []"...
    "Min confinement time for satisfaction of (p*tau_E)_min [s]",...
    "The cost of the powerplant [$]"];
for k = 1:14
    fig = figure;
    y = p1(:,k);
    plot(x, y,'color',rand(1,3));
    title(T(k))
end

