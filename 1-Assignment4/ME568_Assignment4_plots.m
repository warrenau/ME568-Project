%%%%%%%%%%%%%%%%
% Part 1
%%%%%%%%%%%%%%%%

% plots
% plot u, v, w; plot U, V, W; plot the primes
figure(1); clf;
subplot(3,1,1);
pcolor(dat.x, dat.z, u'), shading flat, caxis([-.0025 .0025]),colorbar;
xlabel('x (m)');
ylabel('z (m)');
title("u at t="+time+"s");

subplot(3,1,2);
pcolor(dat.x, dat.z, v'), shading flat, caxis([-.0025 .0025]),colorbar;
xlabel('x (m)');
ylabel('z (m)');
title("v at t="+time+"s");

subplot(3,1,3);
pcolor(dat.x, dat.z, w'), shading flat, caxis([-.0025 .0025]),colorbar;
xlabel('x (m)');
ylabel('z (m)');
title("w at t="+time+"s");

saveas(gcf,'1-plots/Vel_plot_3043.png')

% primes
figure(2); clf;
subplot(3,1,1);
pcolor(dat.x, dat.z, u_prime), shading flat, caxis([-.0025 .0025]),colorbar;
xlabel('x (m)');
ylabel('z (m)');
title("u' at t="+time+"s");

subplot(3,1,2);
pcolor(dat.x, dat.z, v_prime), shading flat, caxis([-.0025 .0025]),colorbar;
xlabel('x (m)');
ylabel('z (m)');
title("v' at t="+time+"s");

subplot(3,1,3);
pcolor(dat.x, dat.z, w_prime), shading flat, caxis([-.0025 .0025]),colorbar;
xlabel('x (m)');
ylabel('z (m)');
title("w' at t="+time+"s");

saveas(gcf,'1-plots/Vel_primes_plot_3043.png')


% U, V, W
figure(3); clf;
subplot(3,1,1);
plot(U,dat.z);
xlabel("U (m/s)");
ylabel("z (m)");
title("U at t="+time+"s");

subplot(3,1,2);
plot(V,dat.z);
xlabel("V (m/s)");
ylabel("z (m)");
title("V at t="+time+"s");

subplot(3,1,3);
plot(W,dat.z);
xlabel("W (m/s)");
ylabel("z (m)");
title("W at t="+time+"s");

saveas(gcf,'1-plots/Vel_avg_plot_3043.png')


% part 1 b

% plot tke; individual and total
figure(4); clf;
subplot(4,1,1);
pcolor(dat.x, dat.z, tke_u), shading flat, colorbar;
xlabel('x (m)');
ylabel('z (m)');
title("u'u' at t="+time+"s");

subplot(4,1,2);
pcolor(dat.x, dat.z, tke_v), shading flat, colorbar;
xlabel('x (m)');
ylabel('z (m)');
title("v'v' at t="+time+"s");

subplot(4,1,3);
pcolor(dat.x, dat.z, tke_w), shading flat, colorbar;
xlabel('x (m)');
ylabel('z (m)');
title("w'w' at t="+time+"s");

subplot(4,1,4);
pcolor(dat.x, dat.z, tke_tot), shading flat, colorbar;
xlabel('x (m)');
ylabel('z (m)');
title("tke at t="+time+"s");

saveas(gcf,'1-plots/tke_plot_3043.png')


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










