% ME 568 Assignment 4
% Austin Warren
% May 2022
clear; clc; clf;

% read in dns data
load dns_data.mat

%%%%%%
% Part 1
%%%%%%
% determine average and fluctuating values; average over x for one z
% coordiate at set time and y values
% we can write a function that will do this for any parameter we choose


% loop set up
% define variables; lets use the second y position (thats what the demo
% used); lets use the middle time position (8); u,v,w should now just be
% functions of x and z
t_ind = 8; %3 for part c of 1
y_ind = 2;

dat = dns_data(t_ind);
time = dat.time;

u = squeeze(dat.u(:,y_ind,:));
v = squeeze(dat.v(:,y_ind,:));
w = squeeze(dat.w(:,y_ind,:));

% loop
numz = dat.nz;
numx = dat.nx;

% u
[U, u_prime] = ReynoldsLoop(u,numz,numx);

% v
[V, v_prime] = ReynoldsLoop(v,numz,numx);

% w
[W, w_prime] = ReynoldsLoop(w,numz,numx);

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



%
% part 1 b
% 
% tke components

tke_u = 0.5*u_prime.^2;
tke_v = 0.5*v_prime.^2;
tke_w = 0.5*w_prime.^2;

tke_tot = (tke_u + tke_v + tke_w);

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

% Reynolds Stress
uu = u_prime.^2;
uv = u_prime.*v_prime;
uw = u_prime.*w_prime;
vv = v_prime.^2;
ww = w_prime.^2;
vw = v_prime.*w_prime;

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



%
% part 1 c
%



%%%%%%
% Part 2
%%%%%%

% part a
% P = - u_i u_j S_ij
% epsilon = 2 nu s_ij s_ij
% S_ij = 1/2 (dU_i/dx_j + dU_j/dx_i)
hx = dat.dx;
hz = dat.dz;

[dUdz] = gradient(U,hz);
[dWdz] = gradient(W,hz);
%[dwdx] = gradient(W,hx,hz);
S_ij = 0.5*(dUdz + dWdz);

[dudz, dudx] = gradient(u_prime, hz, hx);
[dwdz, dwdx] = gradient(w_prime, hz, hx);
[dvdz, dvdx] = gradient(v_prime, hz, hx);
s_ij = 0.5*(dudz + dwdx + dvdx + dvdz);

production = - uw.*S_ij;
dissipation = 2 * dat.nu * (s_ij.*s_ij);

P_horiz = zeros(numz);
ep_horiz = zeros(numz);
for i = 1:numz
   P_horiz(i) = mean(production(i,:)); 
   ep_horiz(i) = mean(dissipation(i,:));
end

P_avg = mean(P_horiz);
ep_avg = mean(ep_horiz);

% plots
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





%
% Functions
%
% Reynolds decomp function
function [U, u_prime] = ReynoldsDecomp(u)
    % u should come in as a vector (one z value, all x values)
    % external loop should loop over z values as needed
    U = mean(u);
    u_prime = u - U;
end
%

% make another function to do the looping
function [U, u_prime] = ReynoldsLoop(u,numz,numx)
    U = zeros(numz,1);
    u_prime = zeros(numz,numx);
    for i = 1:numz
        [U(i), u_prime(i,:)] = ReynoldsDecomp(u(:,i));
    end
end
%


% function for S_ij and s_ij
% function [S_ij] = sij(U, W, hx, hz)
%     [dudx, dudz] = gradient(U,hx,hz); % this might be backwards? x might be the second dimension?
%     [dwdx, dwdz] = gradient(W,hx,hz); % but bc hx=hz and S_ij uses both, it should be fine? -- idk
%     S_ij = 0.5*(dudz + dwdx);
% 
% end
%







