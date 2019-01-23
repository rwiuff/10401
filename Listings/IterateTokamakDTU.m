close all
clear all

titl = ["Desired, output, power, [MW]", "Maximum, wall, load, [MW, m^-2]", ...
    "Magnetic, field, at, the, edge, of, the, coil, [T]", ...
    "Tensile, strenght, of, the, magnetic, field, coils, [atm]"];
foldertitl = ["PE", "PW", "Bmax", "sigmamax"];

for d = 1:4
    l = d; % Choose iteration: 1 = P_E, 2 = P_W, 3 = B_max, 4 = sigma_max
    p0 = [];
    x = [];
    mkdir('../MatlabFigures', foldertitl(l))
    printout = (sprintf('Created folder %s', foldertitl(l)));
    disp(printout)
    if l == 1
        x = linspace(1, 3000, 1000);
        for i = 1:length(x)
            [b, c, a, R_0, A, A_p, V_p, P_dens, p, ...
                n, B_0, beta, tau_E_min, C_per_watt] = ...
                tokamakDTU_asign_1(0.01, 1, 2, x(i), 4, 13, 3000, 0.4);
            q = [b, c, a, R_0, A, A_p, V_p, P_dens, p, ...
                n, B_0, beta, tau_E_min, C_per_watt];
            p0 = cat(1, p0, q);
        end
    elseif l == 2
        x = linspace(1, 10, 1000);
        for i = 1:length(x)
            [b, c, a, R_0, A, A_p, V_p, P_dens, p, ...
                n, B_0, beta, tau_E_min, C_per_watt] = ...
                tokamakDTU_asign_1(0.01, 1, 2, 1000, x(i), 13, 3000, 0.4);
            q = [b, c, a, R_0, A, A_p, V_p, P_dens, p, ...
                n, B_0, beta, tau_E_min, C_per_watt];
            p0 = cat(1, p0, q);
        end
    elseif l == 3
        x = linspace(10, 20, 1000);
        for i = 1:length(x)
            [b, c, a, R_0, A, A_p, V_p, P_dens, p, ...
                n, B_0, beta, tau_E_min, C_per_watt] = ...
                tokamakDTU_asign_1(0.01, 1, 2, 1000, 4, x(i), 3000, 0.4);
            q = [b, c, a, R_0, A, A_p, V_p, P_dens, p, ...
                n, B_0, beta, tau_E_min, C_per_watt];
            p0 = cat(1, p0, q);
        end
    elseif l == 4
        x = linspace(2000, 5000, 1000);
        for i = 1:length(x)
            [b, c, a, R_0, A, A_p, V_p, P_dens, p, ...
                n, B_0, beta, tau_E_min, C_per_watt] = ...
                tokamakDTU_asign_1(0.01, 1, 2, 1000, 4, 13, x(i), 0.4);
            q = [b, c, a, R_0, A, A_p, V_p, P_dens, p, ...
                n, B_0, beta, tau_E_min, C_per_watt];
            p0 = cat(1, p0, q);
        end
    end
    printout = (sprintf('Iterated tokamakDTU_asign_1.m %d of 4', d));
    disp(printout)
    b = p0(:, 1);
    c = p0(:, 2);
    a = p0(:, 3);
    R_0 = p0(:, 4);
    A = p0(:, 5);
    A_p = p0(:, 6);
    V_p = p0(:, 7);
    P_dens = p0(:, 8);
    p = p0(:, 9);
    n = p0(:, 10);
    B_0 = p0(:, 11);
    beta = p0(:, 12);
    beta = beta * 100;
    tau_E_min = p0(:, 13);
    C_per_watt = p0(:, 14);
    
    f1 = figure;
    ptau = plot(x, tau_E_min, 'color', 'r', 'LineWidth', 1);
    title('Min confinement time for satisfaction of (p tau_E)_{min}')
    xlabel(titl(l));
    ylabel('[s]');
    ytickformat('%.2f');
    grid on
    grid minor
    figcount = 1;
    printout = (sprintf('Created figure %d of 8', figcount));
    disp(printout)
    f2 = figure;
    pb = plot(x, b, 'color', 'r', 'LineWidth', 1);
    hold on
    pc = plot(x, c, '--', 'color', 'b', 'LineWidth', 1);
    pa = plot(x, a, '-.', 'color', 'm', 'LineWidth', 1);
    pcpw = plot(x, C_per_watt, 'color', 'g', 'LineWidth', 1);
    hold off
    xlabel(titl(l));
    ytickformat('%.2f');
    lgd = legend([pb, pc, pa, pcpw], ...
        {"Blanket - shield, thickness, [m]", "Magnet, coil, thickness, [m]", ...
        "Minor, radius, [m]", ...
        "Cost, per, watt, [$]"});% Figure legend
    legend('boxoff')
    if l == 3
        legend('Location', 'north')
    else
        legend('Location', 'west')
    end
    grid on
    grid minor
    figcount = 2;
    printout = (sprintf('Created figure %d of 8', figcount));
    disp(printout)
    f3 = figure;
    hold on
    yyaxis left
    pA = plot(x, A, 'LineWidth', 1);
    ylabel('Aspect ratio')
    yyaxis right
    pr_0 = plot(x, R_0, '--', 'LineWidth', 2);
    ylabel('Major radius [m]')
    hold off
    xlabel(titl(l));
    ytickformat('%.2f');
    lgd = legend([pr_0, pA], ...
        {"Major, radius", "Aspect, ratio, []"});% Figure legend
    legend('boxoff')
    legend('Location', 'east')
    grid on
    grid minor
    figcount = 3;
    printout = (sprintf('Created figure %d of 8', figcount));
    disp(printout)
    f4 = figure;
    pA_p = plot(x, A_p, 'LineWidth', 1);
    hold on
    pV_p = plot(x, V_p, '--', 'LineWidth', 2);
    hold off
    xlabel(titl(l));
    ytickformat('%.2f');
    lgd = legend([pA_p, pV_p], ...
        {"Plasma, surface, [m^2]", "Plasma, volume, [m^3]"});% Figure legend
    legend('boxoff')
    legend('Location', 'east')
    grid on
    grid minor
    figcount = 4;
    printout = (sprintf('Created figure %d of 8', figcount));
    disp(printout)
    f5 = figure;
    pP_dens = plot(x, P_dens, 'color', 'r', 'LineWidth', 1);
    title('Power density')
    xlabel(titl(l));
    ylabel('[W m^-1]');
    ytickformat('%.2f');
    grid on
    grid minor
    figcount = 5;
    printout = (sprintf('Created figure %d of 8', figcount));
    disp(printout)
    f6 = figure;
    yyaxis left
    pn = plot(x, n, 'LineWidth', 1);
    ylabel('Particle density [m^-3]')
    ytickformat('%.2f');
    yyaxis right
    p = plot(x, p, '--', 'LineWidth', 2);
    ylabel('Plasma pressure [Pa]')
    ytickformat('%.2f');
    xlabel(titl(l));
    lgd = legend([pn, p], ...
        {"Particle, density", "Plasma, pressure, [Pa]"});% Figure legend
    legend('boxoff')
    legend('Location', 'north')
    grid on
    grid minor
    figcount = 6;
    printout = (sprintf('Created figure %d of 8', figcount));
    disp(printout)
    if l == 1
        f7 = figure;
        subplot(2, 1, 1);
        pB_0 = plot(x, B_0, 'color', 'b', 'LineWidth', 1);
        title('Magnetic field at magnetic axis')
        ylabel('[T]');
        ytickformat('%.2f');
        grid on
        grid minor
        subplot(2, 1, 2);
        B_0(B_0 <= -10) = NaN;
        pB_0 = plot(x, B_0, 'color', 'b', 'LineWidth', 1);
        xlabel(titl(l));
        ylabel('[T]');
        ytickformat('%.2f');
        grid on
        grid minor
    else
        f7 = figure;
        pB_0 = plot(x, B_0, 'color', 'b', 'LineWidth', 1);
        title('Magnetic field at magnetic axis')
        xlabel(titl(l));
        ylabel('[T]');
        ytickformat('%.2f');
        grid on
        grid minor
    end
    figcount = 7;
    printout = (sprintf('Created figure %d of 8', figcount));
    disp(printout)
    if l == 1 || l == 2
        f8 = figure;
        subplot(2, 1, 1);
        pbeta = plot(x, beta, 'color', 'r', 'LineWidth', 1);
        title('Normalised plasma preasure in the centre')
        ylabel('[%]');
        ytickformat('%.2f');
        grid on
        grid minor
        subplot(2, 1, 2);
        beta(beta >= 100) = NaN;
        pbeta = plot(x, beta, 'color', 'r', 'LineWidth', 1);
        xlabel(titl(l));
        ylabel('[%]');
        ytickformat('%.2f');
        grid on
        grid minor
    else
        f8 = figure;
        pbeta = plot(x, beta, 'color', 'r', 'LineWidth', 1);
        title('Normalised plasma preasure in the centre')
        xlabel(titl(l));
        ylabel('[%]');
        ytickformat('%.2f');
        grid on
        grid minor
    end
    figcount = 8;
    printout = (sprintf('Created figure %d of 8', figcount));
    disp(printout)
    % lgd = legend([pb, pc, pa, pr_0, pA, pA_p, pV_p,...
    %     pP_dens, p, pn, pB_0, pbeta, ptau, pcpw], ...
    %     {"Blanket-shield thickness [m]", "Magnet coil thickness [m]"...
    %     "Minor radius [m]", "Major radius [m]", "Aspect ratio []"...
    %     "Plasma surface [m^2]", "Plasma volume [m^3]", "Power density [W m^-1]"...
    %     "Plasma pressure [Pa]", "Particle density [m^-3]",...
    %     "Magnetic field at magnetic axis [T]", "Plasma beta in the centre []"...
    %     "Min confinement time for satisfaction of (p tau_E)_min [s]",...
    %     "The cost of the powerplant [$]"});% Figure legend
    
    
    fignames = ["f1", "f2", "f3", "f4", "f5", "f6", "f7", "f8"];
    figfigs = [f1, f2, f3, f4, f5, f6, f7, f8];
    
    
    for i = 1:8
        foldername = sprintf('../MatlabFigures/%s', foldertitl(l));
        epsfilename = sprintf('%s.eps', fignames(i));
        fullfilename = fullfile(foldername, epsfilename);
        saveas(figfigs(i), fullfilename, 'epsc')
        printout = (sprintf('Saved figure %d of 8', i));
        disp(printout)
    end
    
    close all
    printout = (sprintf('Completed iteration %d of 4', d));
    disp(printout)
end