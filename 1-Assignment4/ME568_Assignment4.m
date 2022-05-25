% ME 568 Assignment 4
% Austin Warren
% May 2022
clear; clc;

% read in dns data
load dns_data.mat

y_ind = 2; % only using the second y plane
P_avg = zeros(length(dns_data),1);
ep_avg = zeros(length(dns_data),1);
J_b_avg = zeros(length(dns_data),1);

for k=1:length(dns_data)
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
    %t_ind = 8; %3 for part c of 1
    

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


    %
    % part 1 b
    % 
    % tke components

    tke_u = 0.5*u_prime.^2;
    tke_v = 0.5*v_prime.^2;
    tke_w = 0.5*w_prime.^2;

    tke_tot = (tke_u + tke_v + tke_w);


    % Reynolds Stress
    uu = u_prime.^2;
    uv = u_prime.*v_prime;
    uw = u_prime.*w_prime;
    vv = v_prime.^2;
    ww = w_prime.^2;
    vw = v_prime.*w_prime;



    %
    % part 1 c
    %
    % lets print out values to a file for the timesteps we are interested
    % in.
    if k==3
        writematrix(u,'2-tabs/u_1354');
        writematrix(v,'2-tabs/v_1354');
        writematrix(w,'2-tabs/w_1354');
        writematrix(U,'2-tabs/Uavg_1354');
        writematrix(V,'2-tabs/Vavg_1354');
        writematrix(W,'2-tabs/Wavg_1354');
        writematrix(u_prime,'2-tabs/uprime_1354');
        writematrix(v_prime,'2-tabs/vprime_1354');
        writematrix(w_prime,'2-tabs/wprime_1354');
        writematrix(tke_u,'2-tabs/utke_1354');
        writematrix(tke_v,'2-tabs/vtke_1354');
        writematrix(tke_w,'2-tabs/wtke_1354');
        writematrix(tke_tot,'2-tabs/tottke_1354');
        writematrix(uu,'2-tabs/uu_1354');
        writematrix(uv,'2-tabs/uv_1354');
        writematrix(uw,'2-tabs/uw_1354');
        writematrix(vv,'2-tabs/vv_1354');
        writematrix(vw,'2-tabs/vw_1354');
        writematrix(ww,'2-tabs/ww_1354');
        
    elseif k==8
        writematrix(u,'2-tabs/u_3043');
        writematrix(v,'2-tabs/v_3043');
        writematrix(w,'2-tabs/w_3043');
        writematrix(U,'2-tabs/Uavg_3043');
        writematrix(V,'2-tabs/Vavg_3043');
        writematrix(W,'2-tabs/Wavg_3043');
        writematrix(u_prime,'2-tabs/uprime_3043');
        writematrix(v_prime,'2-tabs/vprime_3043');
        writematrix(w_prime,'2-tabs/wprime_3043');
        writematrix(tke_u,'2-tabs/utke_3043');
        writematrix(tke_v,'2-tabs/vtke_3043');
        writematrix(tke_w,'2-tabs/wtke_3043');
        writematrix(tke_tot,'2-tabs/tottke_3043');
        writematrix(uu,'2-tabs/uu_3043');
        writematrix(uv,'2-tabs/uv_3043');
        writematrix(uw,'2-tabs/uw_3043');
        writematrix(vv,'2-tabs/vv_3043');
        writematrix(vw,'2-tabs/vw_3043');
        writematrix(ww,'2-tabs/ww_3043');
        
    else
        
    end



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
    S_ij = (dUdz);

    [dudz, dudx] = gradient(u_prime, hz, hx);
    [dwdz, dwdx] = gradient(w_prime, hz, hx);
    [dvdz, dvdx] = gradient(v_prime, hz, hx);


    production = - uw.*S_ij;
    dissipation = zeros(numz,numx);
    for m = 1:numz
        for n = 1:numx
            dissipation(m,n) = 2 * dat.nu * (dudx(m,n).^2 + dudz(m,n).^2 + dvdx(m,n).^2 + dvdz(m,n).^2 + dwdx(m,n).^2 + dwdz(m,n).^2);
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
    
    if k == 8
        writematrix(production,'2-tabs/P_3043');
        writematrix(P_horiz,'2-tabs/P_horiz');
        writematrix(dissipation,'2-tabs/ep_3043');
        writematrix(ep_horiz,'2-tabs/ep_horiz');
    end



    % part 2 b



    %%%%%%
    % Part 3
    %%%%%%

    % part 3 a
    % time integrals of spatially integrated buoyancy production J_b
    % J_b = -g/rho_0 w'rho'

    g = dat.g;
    rho0 = dat.rho0;
    rho = squeeze(dat.density(:,y_ind,:))';
    J_b = -g/rho0 * w_prime.*rho;
    J_b_horiz = zeros(numz,1);
    for j=1:numz
       J_b_horiz(j) = mean(J_b(j,:)); 
    end
    J_b_avg(k) = mean(J_b_horiz);
    Gamma_avg = J_b_avg ./ ep_avg;



end

% write out produciton, dissipation, and buoyancy at each time
writematrix(P_avg,'2-tabs/production');
writematrix(ep_avg,'2-tabs/dissipation');
writematrix(J_b_avg,'2-tabs/buoyancy');
writematrix(Gamma_avg,'2-tabs/mixing');


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










