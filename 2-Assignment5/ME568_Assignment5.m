% ME 568 Assignment 5
% Austin Warren
% May 2022
clear; clc; %figure(1);clf;


% read in dns data
load dns_data.mat

y_ind = 2; % only using the second y plane
P_avg = zeros(length(dns_data),1);
ep_avg = zeros(length(dns_data),1);
P_sum = zeros(length(dns_data),1);
ep_sum = zeros(length(dns_data),1);
J_b_avg = zeros(length(dns_data),1);
tke_sum = zeros(1,length(dns_data));
U_time = zeros(length(dns_data), dns_data(1).nz);

char_ell = zeros(dns_data(1).nz,length(dns_data));
lambda = zeros(dns_data(1).nz,length(dns_data));
eta_horiz = zeros(dns_data(1).nz,length(dns_data));

spectraz = 10;
nfft = 256;
power = zeros(length(dns_data),nfft/2);
freq = zeros(length(dns_data),nfft/2);
for k=1:length(dns_data)


    dat = dns_data(k);
    time = dat.time;
    
    % pull out velocities of interest
    u = squeeze(dat.u(:,y_ind,:));
    v = squeeze(dat.v(:,y_ind,:));
    w = squeeze(dat.w(:,y_ind,:));

    %
    numz = dat.nz;
    numx = dat.nx;
    
    % Reynolds decomp -- MATLAB requires functions at the end of the file
    % u
    [U, u_prime] = ReynoldsLoop(u,numz,numx);
    U_time(k,:) = U;
    % v
    [V, v_prime] = ReynoldsLoop(v,numz,numx);

    % w
    [W, w_prime] = ReynoldsLoop(w,numz,numx);
    
    
    % tke components

    tke_u = 0.5*u_prime.^2;
    tke_v = 0.5*v_prime.^2;
    tke_w = 0.5*w_prime.^2;

    tke_tot = (tke_u + tke_v + tke_w);
    tke_sum(k) = sum(tke_tot,'all');
    
    %
    % Reynolds Stress
    uu = u_prime.^2;
    uv = u_prime.*v_prime;
    uw = u_prime.*w_prime;
    vv = v_prime.^2;
    ww = w_prime.^2;
    vw = v_prime.*w_prime;

    %
    hx = dat.dx;
    hz = dat.dz;

    [dUdz] = gradient(U,hz);
    S_ij = (dUdz);

    [dudz, dudx] = gradient(u_prime, hz, hx);
    [dwdz, dwdx] = gradient(w_prime, hz, hx);
    [dvdz, dvdx] = gradient(v_prime, hz, hx);


    production = - uw.*S_ij;
    dissipation = zeros(numz,numx);
    for m = 1:numz
        for n = 1:numx
            sij = [0.5*(dudx(m,n) + dudx(m,n)), 0.5*(dudz(m,n) + dwdx(m,n));
                   0.5*(dwdx(m,n) + dudz(m,n)), 0.5*(dwdz(m,n) + dwdz(m,n))]; % need to add dy terms, but I dont want to add a third dimension to the code right now
            dissipation(m,n) = 2 * dat.nu * sum(sij.^2,'all');
        end
    end

    P_horiz = zeros(numz,1);
    ep_horiz = zeros(numz,1);
    for i = 1:numz
       P_horiz(i,1) = mean(production(i,:));
       ep_horiz(i,1) = mean(dissipation(i,:));
    end

    P_avg(k) = mean(P_horiz);
    ep_avg(k) = mean(ep_horiz);
    
    P_sum(k) = sum(production,'all');
    ep_sum(k) = sum(dissipation,'all');
    
    % kolmogorov
    eta = abs(dat.nu ./ u_prime); 
    % correlation for integral length scale (and Taylor length scale?)
    
    maxlag = numx-2;
    for i=1:numz
        [Rxy, rhoxy, s2x, s2y, mux, muy, lag, Nk] = xcovar(u_prime(i,:),u_prime(i,:),maxlag);
        char_ell(i,k) = trapz(dat.dx, rhoxy(numx+1:end));
        drhoxy = diff(rhoxy(numx+1:end))/dat.dx;
        ddrhoxy = diff(drhoxy)/dat.dx;
        lambda(i,k) = real(sqrt(-2/ddrhoxy(1)));
        eta_horiz(i,k) = mean(eta(i,:));
        
        if k==8 && (i==ceil(numz/2) || i==ceil(numz/4) || i==ceil(3*numz/4) )
            %figure(1);
            %plot(lag, rhoxy);
            %hold on
        end
    end
    
    % spectra stuff -- average over first 10 or so z slices?
    p = zeros(nfft/2,spectraz);
    f = zeros(nfft/2,spectraz);
    for j=1:spectraz
        [p(:,j),f(:,j)] = fast_psd(u_prime(j,:),nfft,1);
    end
    
    for h=1:nfft/2
        power(k,h) = mean(p(h,:));
        freq(k,h) = mean(f(h,:));
    end
    
    

    
    
    
    
    
    
end


%figure(1);
%xlabel('lag')
%ylabel('\rho')
%legend('i=1/2','i=1/4','i=3/4')

% plot psd
figure(2); clf;
for i=1:length(dns_data)
    loglog(freq(i,:),power(i,:))
    hold on
end
xlabel('Frequency')
ylabel('\phi')
legend('671.8','1010.5','1354.1','1687.6','2031.0','2369.9','2704.6','3043.9','3383.3','3722.7','4067.0','4406.3','4749.6','5291.5','5827.2','6505.0','location','southwest','NumColumns',2)
saveas(gcf,'1-plots/psd_plot.png')

% characteristic velocity
char_vel = sqrt(tke_sum);
% characteristic length
char_length = max(abs(char_ell));
% characteristic time scale
char_time = char_length./char_vel;
% dissipation from characteristic velocity and length
char_diss = char_vel.^3 ./ char_length;
% turbulent viscosity based on characteristic velocity and length
char_nuT = char_vel .* char_length;
nuT = mean(char_nuT);
U_nuT = zeros(length(dns_data), dns_data(1).nz);
U_nuT(1,:) = U_time(1,:);
time_int = zeros(length(dns_data),1);
dUdz_time = zeros(length(dns_data)-1, dns_data(1).nz);

for i=length(dns_data)-1
   dUdz_time(i,:) = gradient(U_time(i,:),dns_data(1).dx);
   U_nuT(i+1,:) =  U_nuT(i) + (time_int(i+1)-time_int(i))*nuT*dUdz_time(i,:);
end

% plot U(z) for nuT verification
figure(1); clf;
plot(U_time(1,:), dns_data(1).z,'-k','linewidth',2)
hold on
plot(U_time(end,:), dns_data(end).z,'--k','linewidth',2)
plot(U_nuT(end,:), dns_data(end).z,'--b','linewidth',2)
xlabel('Velocity (m/s)')
ylabel('z (m)')
legend('U_{init}','U_{final}','U_{final,\nu}')



% integrate dissipation and production

for i=1:length(dns_data)
    time_int(i) = dns_data(i).time;
end
prod_int = trapz(time_int, P_avg);
diss_int = trapz(time_int, ep_avg);

P_sum_tot = sum(P_sum,'all');
ep_sum_tot = sum(ep_sum,'all');

% compare length scales
% take the median at each time of the Taylor and Komogorov scales
lambda_med = zeros(length(dns_data),1);
eta_med = zeros(length(dns_data),1);
for i=1:length(dns_data)
    lambda_med(i) = median(lambda(:,i));
    eta_med(i) = median(eta_horiz(:,i));
end

% Kolmogorov calcs with dissipation
eta_diss = (dns_data(1).nu^3 ./ ep_sum).^(1/4);
eta_time_diss = (dns_data(1).nu ./ ep_sum).^(1/2);
eta_vel_diss = (dns_data(1).nu .* ep_sum).^(1/4);


% plot length scales
figure(3); clf;
semilogy(time_int,char_length,'-k','linewidth',2)
hold on
semilogy(time_int,lambda_med,'-b','linewidth',2)
semilogy(time_int,eta_med,'-r','linewidth',2)
semilogy(time_int,eta_diss,'--r','linewidth',2)
xlabel('Time (s)')
ylabel('Length (m)')
legend('$\ell$','$\lambda$','$\eta$','$\eta_{\epsilon}$','interpreter','latex','location','northoutside')
%ylim([1e-3 1e1]);
saveas(gcf,'1-plots/char_length_comp_plot.png')

% plot time scales
figure(4); clf;
semilogy(time_int,char_time,'-k','linewidth',2)
hold on
semilogy(time_int,eta_time_diss,'-r','linewidth',2)
xlabel('Time (s)')
ylabel('Time scale (s)')
legend('$t$','$\tau_{\eta}$','interpreter','latex','location','east')
saveas(gcf, '1-plots/char_time_plot.png')

% plot eddy viscosity
figure(5); clf;
plot(time_int, char_nuT, '-k','linewidth',2)
xlabel('Time (s)')
ylabel('\nu_T (m^2/s)')
saveas(gcf, '1-plots/char_nuT_plot.png')

% plot velocity scales
figure(6); clf;
semilogy(time_int, char_vel,'-k','linewidth',2)
hold on
semilogy(time_int, eta_vel_diss, '-r','linewidth',2)
xlabel('Time (s)')
ylabel('Velocity (m/s)')
legend('$u_{char}$','$u_{\eta}$','interpreter','latex')
saveas(gcf, '1-plots/char_vel_plot.png')






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
%
% auto-correlation function






