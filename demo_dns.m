% demo file to illustrate how to read/interpret dns data.
%
% JDN, April 15, 2022

load dns_data.mat

% dns_data is a 1x16 structure containing many variables:
% 
% dns_data(n) % this returns the nth time-step
% 
% ans = 
% 
%        time: 604.1116
%           x: [384x1 double]
%           y: [3x1 double]
%           z: [201x1 double]
%      p_phys: [20x1 double]
%        p_ic: [20x1 double]
%           u: [384x3x201 double]
%           v: [384x3x201 double]
%           w: [384x3x201 double]
%     density: [384x3x201 double]
%          nx: 384
%          ny: 80
%          nz: 201
%          dx: 0.0109
%          dy: 0.0109
%          dz: 0.0109
%          Lx: 4.1888
%          Ly: 0.8727
%          Lz: 2.1817
%         h_0: 0.3000
%         k_0: 3.0000
%         d_u: 0.0106
%         d_r: 3.0268e-09
%         z_0: 1.0908
%         rsc: 1
%           a: 0.2000
%           b: 0.5000
%           c: 0.1000
%           g: 9.8100
%        rho0: 1027
%          nu: 1.0000e-06
%       kappa: 5.7143e-07
%        d_rr: 0.0031
% 


% to plot an x-z section of density, du/dx, du/dy and du/dz at each time step for
% the middle y-index (y_ind=2)

for a=1:length(dns_data)
    
    dat=dns_data(a);
    
    % first compute gradients:
    
    % NOTE gradient for some reason returns the gradients in a transposed order:
    % [d/dy, d/dx, d/dz]  ****
    [dudy,dudx,dudz]=gradient(dat.u,dat.dx); % this returns the three components of the 
                                             % gradient of the x-component of velocity
                                             % note that this form
                                             % only works because dx=dy=dz
    y_ind=2; % the y-index we will plot
    figure(1)
    subplot(4,1,1)
    pcolor(dat.x,dat.z,squeeze(dat.density(:,y_ind,:))'), shading flat, caxis([-.0025 .0025]),colorbar
    text(0.02,1.1,'density anomaly [kg/m^3]','units','normalized')
    title(['time = ' num2str(dat.time) ' s'])
    % squeeze changes the dimension from [384,1,201] to [384,201]
    % the prime is needed because of matlab's row,column conventions
    
    % note that you need to add dat.rho0 (1027 kg/m^3) to get the real
    % density
    
    clims=[-.05 .05]; % limits for plotting
    subplot(4,1,2)
    pcolor(dat.x,dat.z,squeeze(dudx(:,y_ind,:))'), shading flat, caxis(clims),colorbar
    text(0.02,.9,'du/dx [s^{-1}]','units','normalized')
    subplot(4,1,3)
    pcolor(dat.x,dat.z,squeeze(dudy(:,y_ind,:))'), shading flat, caxis(clims),colorbar
    text(0.02,.9,'du/dy [s^{-1}]','units','normalized')
    subplot(4,1,4)
    pcolor(dat.x,dat.z,squeeze(dudz(:,y_ind,:))'), shading flat, caxis(clims),colorbar
    xlabel('X [m]'),ylabel('Z [m]')
    text(0.02,.9,'du/dz [s^{-1}]','units','normalized')

    colormap(redblue4)
    pause
end