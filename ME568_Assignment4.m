% ME 568 Assignment 4
% Austin Warren
% May 2022

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

