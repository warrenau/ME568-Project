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
t_ind = 8;
y_ind = 2;

dat = dns_data(t_ind);

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


%
% part 1 b
% 
% tke components
tke = zeros(3,numz,numx);
tke(1,:,:) = u_prime.^2;
tke(2,:,:) = v_prime.^2;
tke(3,:,:) = w_prime.^2;

% Reynolds Stress
uv = u_prime.*v_prime;
uw = u_prime.*w_prime;

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

[dudz] = gradient(U,hz); 
%[dwdx] = gradient(W,hx,hz);
S_ij = 0.5*(dudz);

[dudz, dudx] = gradient(u_prime, hz, hx);
[dwdz, dwdx] = gradient(w_prime, hz, hx);
s_ij = 0.5*(dudz + dwdx);

production = - uw.*S_ij;
dissipation = 2 * dat.nu * (s_ij.*s_ij);





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







