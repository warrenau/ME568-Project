% ME 568 Assignment 5
% Austin Warren
% May 2022
clear; clc; clf;


% read in dns data
load dns_data.mat

y_ind = 2; % only using the second y plane
P_avg = zeros(length(dns_data),1);
ep_avg = zeros(length(dns_data),1);
J_b_avg = zeros(length(dns_data),1);
tke_sum = zeros(length(dns_data),1);


char_ell = zeros(dns_data(1).nz,length(dns_data));

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
    
    
    % correlation for integral length scale (and Taylor length scale?)
    
    maxlag = numx;
    for i=1:numz
        [Rxy, rhoxy, s2x, s2y, mux, muy, lag, Nk] = xcovar(u_prime,u_prime,maxlag);
        char_ell(i,k) = trapz(rhoxy,lag);
        
        if k==8 && i==ceil(numz/2)
            hold on
            plot(lag, rhoxy);
            xlabel('lag')
            ylabel('\rho')
        end
    end
    
    

    
    
    
    
    
    
end


% characteristic velocity
char_vel = sqrt(tke_sum);




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






