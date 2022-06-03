% ME 568 Assignment 5
% Austin Warren
% May 2022
clear; clc; clf;


% read in dns data
load dns_data.mat

y_ind = 2; % only using the second y plane
P_avg = zeros(length(dns_data),1);
ep_avg = zeros(length(dns_data),1);
P_sum = zeros(length(dns_data),1);
ep_sum = zeros(length(dns_data),1);
J_b_avg = zeros(length(dns_data),1);
tke_sum = zeros(1,length(dns_data));


char_ell = zeros(dns_data(1).nz,length(dns_data));
lambda = zeros(dns_data(1).nz,length(dns_data));
eta_horiz = zeros(dns_data(1).nz,length(dns_data));

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
            plot(lag, rhoxy);
            hold on
            xlabel('lag')
            ylabel('\rho')
        end
    end
    
    

    
    
    
    
    
    
end

legend('i=1/2','i=1/4','i=3/4')
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


% integrate dissipation and production
time_int = zeros(length(dns_data),1);
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
figure(2); clf;
semilogy(time_int,char_length,'-k','linewidth',2)
hold on
semilogy(time_int,lambda_med,'-b','linewidth',2)
semilogy(time_int,eta_med,'-r','linewidth',2)
xlabel('Time (s)')
ylabel('Length (m)')
legend('$\ell$','$\lambda$','$\eta$','interpreter','latex')
ylim([1e-3 1e1]);
saveas(gcf,'1-plots/char_length_comp_plot.png')





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






