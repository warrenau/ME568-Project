% ME 568 Assignment 4
% Austin Warren
% May 2022

% read in dns data
load dns_data.mat

% determine average and fluctuating values; average over x for one z
% coordiate at set time and y values
% we can write a function that will do this for any parameter we choose
function [U, u_prime] = ReynoldsDecomp(u)
    % u should come in as a vector (one z value, all x values)
    % external loop should loop over z values as needed
    U = mean(u);
    u_prime = u - U;
end


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
U = zeros(numz,1);
u_prime = zeros(numz,numx);

for i = 1:numz
    [U(i),u_prime(i)] = ReynoldsDecomp(u);





end