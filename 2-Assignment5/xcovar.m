function [Rxy, rhoxy, s2x, s2y, mux, muy, lag, Nk] = xcovar(x,y,maxlag)

    lag = -maxlag:1:maxlag;
    Ntot = length(x);

    ctr = 0;

    % calculate negative and zero lag
    for k = -maxlag:1:0
        ctr = ctr + 1; % update counter

        nn1 = 1+abs(k):1:Ntot;
        nn2 = 1:1:Ntot-abs(k);

        % setup temporary time series 
        xx = x(nn1); % x(t)
        yy = y(nn2); % y(t - tau)

        % make sure we are only using data values that aren't missing in
        % either xx or yy

        nn = find(~isnan(xx.*yy));
        Nk(ctr) = length(nn);
        mux(ctr) = mean(xx(nn));
        s2x(ctr) = sum((xx(nn)-mux(ctr)).^2)/(Nk(ctr)-1);
        muy(ctr) = mean(yy(nn));
        s2y(ctr) = sum((yy(nn)-muy(ctr)).^2)/(Nk(ctr)-1);
        Rxy(ctr) = sum((xx(nn)-mux(ctr)).*(yy(nn)-muy(ctr)))/(Nk(ctr)-1);
        rhoxy(ctr) = Rxy(ctr)/sqrt(s2x(ctr)*s2y(ctr));
    end

    % calculate positive lags
    for k = 1:1:maxlag
        ctr = ctr + 1;
        nn1 = 1:1:Ntot-k;
        nn2 = k+1:1:Ntot;

         % setup temporary time series 
        xx = x(nn1); % x(t)
        yy = y(nn2); % y(t - tau)

        % make sure we are only using data values that aren't missing in
        % either xx or yy

        nn = find(~isnan(xx.*yy));
        Nk(ctr) = length(nn);
        mux(ctr) = mean(xx(nn));
        s2x(ctr) = sum((xx(nn)-mux(ctr)).^2)/(Nk(ctr)-1);
        muy(ctr) = mean(yy(nn));
        s2y(ctr) = sum((yy(nn)-muy(ctr)).^2)/(Nk(ctr)-1);
        Rxy(ctr) = sum((xx(nn)-mux(ctr)).*(yy(nn)-muy(ctr)))/(Nk(ctr)-1);
        rhoxy(ctr) = Rxy(ctr)/sqrt(s2x(ctr)*s2y(ctr));
    end
end