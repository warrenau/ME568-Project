% plotting script for ME 568 Assignment 4
clear; clc;
load dns_data.mat

dat = dns_data(3); % pick a set of dns_data to use for dx, dz, etc

%%%%%%%%%%%%%%%%
% Part 1
%%%%%%%%%%%%%%%%

% time = 1354.1 seconds
U_1354 = readmatrix('2-tabs/Uavg_1354');
V_1354 = readmatrix('2-tabs/Vavg_1354');
W_1354 = readmatrix('2-tabs/Wavg_1354');
u_1354 = readmatrix('2-tabs/u_1354');
v_1354 = readmatrix('2-tabs/v_1354');
w_1354 = readmatrix('2-tabs/w_1354');
uprime_1354 = readmatrix('2-tabs/uprime_1354');
vprime_1354 = readmatrix('2-tabs/vprime_1354');
wprime_1354 = readmatrix('2-tabs/wprime_1354');
tke_u_1354 = readmatrix('2-tabs/utke_1354');
tke_v_1354 = readmatrix('2-tabs/vtke_1354');
tke_w_1354 = readmatrix('2-tabs/wtke_1354');
tke_tot_1354 = readmatrix('2-tabs/tottke_1354');
uu_1354 = readmatrix('2-tabs/uu_1354');
uv_1354 = readmatrix('2-tabs/uv_1354');
uw_1354 = readmatrix('2-tabs/uw_1354');
vv_1354 = readmatrix('2-tabs/vv_1354');
vw_1354 = readmatrix('2-tabs/vw_1354');
ww_1354 = readmatrix('2-tabs/ww_1354');

% time = 3043.9 seconds
U_3043 = readmatrix('2-tabs/Uavg_3043');
V_3043 = readmatrix('2-tabs/Vavg_3043');
W_3043 = readmatrix('2-tabs/Wavg_3043');
u_3043 = readmatrix('2-tabs/u_3043');
v_3043 = readmatrix('2-tabs/v_3043');
w_3043 = readmatrix('2-tabs/w_3043');
uprime_3043 = readmatrix('2-tabs/uprime_3043');
vprime_3043 = readmatrix('2-tabs/vprime_3043');
wprime_3043 = readmatrix('2-tabs/wprime_3043');
tke_u_3043 = readmatrix('2-tabs/utke_3043');
tke_v_3043 = readmatrix('2-tabs/vtke_3043');
tke_w_3043 = readmatrix('2-tabs/wtke_3043');
tke_tot_3043 = readmatrix('2-tabs/tottke_3043');
uu_3043 = readmatrix('2-tabs/uu_3043');
uv_3043 = readmatrix('2-tabs/uv_3043');
uw_3043 = readmatrix('2-tabs/uw_3043');
vv_3043 = readmatrix('2-tabs/vv_3043');
vw_3043 = readmatrix('2-tabs/vw_3043');
ww_3043 = readmatrix('2-tabs/ww_3043');


% plots
% plot u, v, w; plot U, V, W; plot the primes

u_names = {u_3043, u_1354, v_3043, v_1354, w_3043, w_1354};
times = ["3043.9","1354.1","3043.9","1354.1","3043.9","1354.1"];
figure(1); clf;
for k = 1:length(u_names)
    subplot(3,2,k);
    pcolor(dat.x, dat.z, u_names{k}'), shading flat, caxis([-.0025 .0025]),colorbar;
    xlabel('x (m)');
    ylabel('z (m)');
    title("u at t="+times(k)+"s");
end
saveas(gcf,'1-plots/Vel_plot_3043_1354.png')

% primes
prime_names = {uprime_3043,uprime_1354,vprime_3043,vprime_1354,wprime_3043,wprime_1354};
figure(2); clf;
for k=1:length(prime_names)
    subplot(3,2,k);
    pcolor(dat.x, dat.z, prime_names{k}), shading flat, caxis([-.0025 .0025]),colorbar;
    xlabel('x (m)');
    ylabel('z (m)');
    title("u' at t="+times(k)+"s");
end
saveas(gcf,'1-plots/Vel_primes_plot_3043_1354.png')


% U, V, W
U_names = {U_3043, U_1354, V_3043, V_1354, W_3043, W_1354};
figure(3); clf;
for k=1:length(U_names)
    subplot(3,2,k);
    plot(U_names{k},dat.z);
    xlabel("U (m/s)");
    ylabel("z (m)");
    title("U at t="+times(k)+"s");
end
saveas(gcf,'1-plots/Vel_avg_plot_3043_1354.png')


% part 1 b

% plot tke; individual and total
tke_names = {tke_u_3043, tke_u_1354, tke_v_3043, tke_v_1354, tke_w_3043, tke_w_1354, tke_tot_3043, tke_tot_1354};
tke_times = ["3043.9","1354.1","3043.9","1354.1","3043.9","1354.1","3043.9","1354.1"];
figure(4); clf;
for k=1:length(tke_names)
    subplot(4,2,k);
    pcolor(dat.x, dat.z, tke_names{k}), shading flat, colorbar;
    xlabel('x (m)');
    ylabel('z (m)');
    title("u'u' at t="+tke_times(k)+"s");
end
saveas(gcf,'1-plots/tke_plot_3043_1354.png')


uu_names_3043 = {uu_3043,uv_3043,uw_3043,vv_3043,vw_3043,ww_3043};
figure(6); clf;
for k=1:6
    subplot(3,2,k);
    pcolor(dat.x, dat.z, uu_names_3043{k}), shading flat, colorbar;
    xlabel("x (m)");
    ylabel("z (m)");
    title("u'u' at t=3043.9s");
end
saveas(gcf,'1-plots/ReynoldsStress_plot_3043.png')

uu_names_1354 = {uu_1354,uv_1354,uw_1354,vv_1354,vw_1354,ww_1354};
figure(7); clf;
for k=1:6
    subplot(3,2,k);
    pcolor(dat.x, dat.z, uu_names_1354{k}), shading flat, colorbar;
    xlabel("x (m)");
    ylabel("z (m)");
    title("u'u' at t=1354.1s");
end
saveas(gcf,'1-plots/ReynoldsStress_plot_1354.png')

% part 1 c
% see above


%%%%%%%%%%%%%%%%
% Part 2
%%%%%%%%%%%%%%%%
P_avg = readmatrix('2-tabs/production');
P_horiz = readmatrix('2-tabs/P_horiz');
P_3043 = readmatrix('2-tabs/P_3043');
ep_avg = readmatrix('2-tabs/dissipation');
ep_horiz = readmatrix('2-tabs/ep_horiz');
ep_3043 = readmatrix('2-tabs/ep_3043');

% plots production and dissipation
figure(8); clf;
subplot(2,1,1);
pcolor(dat.x, dat.z, P_3043), shading flat, colorbar;
xlabel('x (m)');
ylabel('z (m)');
title("P at t=3043.9s");

subplot(2,1,2);
pcolor(dat.x, dat.z, ep_3043), shading flat, colorbar;
xlabel('x (m)');
ylabel('z (m)');
title("\epsilon at t=3043.9s");

saveas(gcf,'1-plots/prod_diss_plot_3043.png')

%
figure(9);clf;
subplot(2,1,1);
plot(P_horiz, dat.z);
xlabel("P");
ylabel("z (m)");
title("Horizontal average P at t=3043.9s")

subplot(2,1,2);
plot(P_horiz, dat.z);
xlabel("\epsilon");
ylabel("z (m)");
title("Horizontal average \epsilon at t=3043.9s")

saveas(gcf,'1-plots/Horizontal_prod_diss_plot_3043.png')

% part 2 b
% see below, plotting with J_b





%%%%%%%%%%%%%%%%
% Part 3
%%%%%%%%%%%%%%%%

J_b_avg = readmatrix('2-tabs/buoyancy');

figure(10); clf;
subplot(3,1,1);
plot(dns_data.time, P_avg);
xlabel("time (s)");
ylabel("P");
title("Volume average P over time")

subplot(3,1,2);
plot(dns_data.time, ep_avg);
xlabel("time (s)");
ylabel("\epsilon");
title("Volume average \epsilon over time");

subplot(3,1,3);
plot(dns_data.time, J_b_avg);
xlabel("time (s)");
ylabel("J_b");
title("Volume average J_b over time");

saveas(gcf,'1-plots/Time_prod_diss_Jb_plot.png')








