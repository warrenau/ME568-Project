% plotting script for ME 568 Assignment 4
load dns_data.m

dat = dns_data(3); % pick a set of dns_data to use for dx, dz, etc

%%%%%%%%%%%%%%%%
% Part 1
%%%%%%%%%%%%%%%%

% time = 1354.1 seconds
U_1354 = readmatrix("2-tabs/U_1354.csv");
V_1354 = readmatrix("2-tabs/V_1354.csv");
W_1354 = readmatrix("2-tabs/W_1354.csv");
uprime_1354 = readmatrix("2-tabs/uprime_1354.csv");
vprime_1354 = readmatrix("2-tabs/vprime_1354.csv");
wprime_1354 = readmatrix("2-tabs/wprime_1354.csv");
tke_u_1354 = readmatrix("2-tabs/utke_1354.csv");
tke_v_1354 = readmatrix("2-tabs/vtke_1354.csv");
tke_w_1354 = readmatrix("2-tabs/wtke_1354.csv");
tke_tot_1354 = readmatrix("2-tabs/tottke_1354.csv");
uu_1354 = readmatrix("2-tabs/uu_1354.csv");
uv_1354 = readmatrix("2-tabs/uv_1354.csv");
uw_1354 = readmatrix("2-tabs/uw_1354.csv");
vv_1354 = readmatrix("2-tabs/vv_1354.csv");
vw_1354 = readmatrix("2-tabs/vw_1354.csv");
ww_1354 = readmatrix("2-tabs/ww_1354.csv");

% time = 3043.9 seconds
U_3043 = readmatrix("2-tabs/U_3043.csv");
V_3043 = readmatrix("2-tabs/V_3043.csv");
W_3043 = readmatrix("2-tabs/W_3043.csv");
uprime_3043 = readmatrix("2-tabs/uprime_3043.csv");
vprime_3043 = readmatrix("2-tabs/vprime_3043.csv");
wprime_3043 = readmatrix("2-tabs/wprime_3043.csv");
tke_u_3043 = readmatrix("2-tabs/utke_3043.csv");
tke_v_3043 = readmatrix("2-tabs/vtke_3043.csv");
tke_w_3043 = readmatrix("2-tabs/wtke_3043.csv");
tke_tot_3043 = readmatrix("2-tabs/tottke_3043.csv");
uu_3043 = readmatrix("2-tabs/uu_3043.csv");
uv_3043 = readmatrix("2-tabs/uv_3043.csv");
uw_3043 = readmatrix("2-tabs/uw_3043.csv");
vv_3043 = readmatrix("2-tabs/vv_3043.csv");
vw_3043 = readmatrix("2-tabs/vw_3043.csv");
ww_3043 = readmatrix("2-tabs/ww_3043.csv");


% plots
% plot u, v, w; plot U, V, W; plot the primes

u_names = [u_3043, u_1354, v_3043, v_1354, w_3043, w_1354];
times = ["3043.9","1354.1","3043.9","1354.1","3043.9","1354.1"];
figure(1); clf;
for k = 1:length(u_names)
    subplot(3,2,k);
    pcolor(dat.x, dat.z, u_names(k)'), shading flat, caxis([-.0025 .0025]),colorbar;
    xlabel('x (m)');
    ylabel('z (m)');
    title("u at t="+times(k)+"s");
end
saveas(gcf,'1-plots/Vel_plot_3043_1354.png')

% primes
prime_names = [uprime_3043,uprime_1354,vprime_3043,vprime_1354,wprime_3043,wprime_1354];
figure(2); clf;
for k=1:length(prime_names)
    subplot(3,2,k);
    pcolor(dat.x, dat.z, prime_names(k)), shading flat, caxis([-.0025 .0025]),colorbar;
    xlabel('x (m)');
    ylabel('z (m)');
    title("u' at t="+times(k)+"s");
end
saveas(gcf,'1-plots/Vel_primes_plot_3043_1354.png')


% U, V, W
U_names = [U_3043, U_1354, V_3043, V_1354, W_3043, W_1354];
figure(3); clf;
for k=1:length(U_names)
    subplot(3,2,k);
    plot(U_names(k),dat.z);
    xlabel("U (m/s)");
    ylabel("z (m)");
    title("U at t="+times(k)+"s");
end
saveas(gcf,'1-plots/Vel_avg_plot_3043_1354.png')


% part 1 b

% plot tke; individual and total
tke_names = [tke_u_3043, tke_u_1354, tke_v_3043, tke_v_1354, tke_w_3043, tke_w_1354];
tke_times = ["3043.9","1354.1","3043.9","1354.1","3043.9","1354.1","3043.9","1354.1"];
figure(4); clf;
for k=1:length(tke_names)
    subplot(4,2,k);
    pcolor(dat.x, dat.z, tke_names(k)), shading flat, colorbar;
    xlabel('x (m)');
    ylabel('z (m)');
    title("u'u' at t="+tke_times(k)+"s");
end
saveas(gcf,'1-plots/tke_plot_3043_1354.png')


figure(6); clf;
subplot(3,2,1);
pcolor(dat.x, dat.z, uu), shading flat, colorbar;
xlabel("x (m)");
ylabel("z (m)");
title("u'u' at t="+time+"s");

subplot(3,2,2);
pcolor(dat.x, dat.z, uv), shading flat, colorbar;
xlabel("x (m)");
ylabel("z (m)");
title("u'v' at t="+time+"s");

subplot(3,2,3);
pcolor(dat.x, dat.z, uw), shading flat, colorbar;
xlabel("x (m)");
ylabel("z (m)");
title("u'w' at t="+time+"s");

subplot(3,2,4);
pcolor(dat.x, dat.z, vv), shading flat, colorbar;
xlabel("x (m)");
ylabel("z (m)");
title("v'v' at t="+time+"s");

subplot(3,2,5);
pcolor(dat.x, dat.z, vw), shading flat, colorbar;
xlabel("x (m)");
ylabel("z (m)");
title("v'w' at t="+time+"s");

subplot(3,2,6);
pcolor(dat.x, dat.z, ww), shading flat, colorbar;
xlabel("x (m)");
ylabel("z (m)");
title("w'w' at t="+time+"s");

saveas(gcf,'1-plots/ReynoldsStress_plot_3043.png')


% part 1 c



%%%%%%%%%%%%%%%%
% Part 2
%%%%%%%%%%%%%%%%

% plots production and dissipation
figure(7); clf;
subplot(2,1,1);
pcolor(dat.x, dat.z, production), shading flat, colorbar;
xlabel('x (m)');
ylabel('z (m)');
title("P at t="+time+"s");

subplot(2,1,2);
pcolor(dat.x, dat.z, dissipation), shading flat, colorbar;
xlabel('x (m)');
ylabel('z (m)');
title("\epsilon at t="+time+"s");

saveas(gcf,'1-plots/prod_diss_plot_3043.png')

% part 2 b




%%%%%%%%%%%%%%%%
% Part 3
%%%%%%%%%%%%%%%%










