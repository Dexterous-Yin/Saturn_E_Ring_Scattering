% calc collision time and frequency
function [coll_num,coll_length] = calc_coll(L,PA,params,ringparams)
% PA: 1xn in deg; L: 1x1
coll_num = zeros(length(ringparams.size),length(PA))+nan;
coll_length = zeros(1,length(PA))+nan;
L_now = L;
lc_eq = asin(sqrt(1/sqrt(4*L_now^6-3*L_now^5)));
for i=1:length(PA)
    pa_now = PA(i);
    if pa_now>=0 %lc_eq/pi*180
        if abs(pa_now-90)<1e-3
            X_now = [L_now;0;0];
            track_length = 4.*L.*params.Rs.*(1.3802-0.3198.*(sind(pa_now)+sind(pa_now).^(1/2)));
            ringn = ringparams.ringn.*get_ring_n(X_now);
            ringnum = ringparams.cross.*track_length.*ringn;
            coll_length(i) = track_length;
            coll_num(:,i) = ringnum;
        else
            F = @(x) params.Bs*1e9/L_now^3*sqrt(1+3*x^2)/(1-x^2)^3-params.Bs*1e9/L_now^3/sind(pa_now)^2;
            sm = abs(fzero(F,[0,0.999]));
            deltas = sm/100000;
            sdots = -sm+deltas/2:deltas:sm-deltas/2;
            Bdots = sqrt(1+3.*sdots.^2)./(1-sdots.^2).^3;
            PAdots = asin(sqrt(Bdots.*sind(pa_now)^2));
            Lengthdots = calc_field_length(L_now,sdots+deltas/2)-calc_field_length(L_now,sdots-deltas/2);
            totalLengthdots = Lengthdots*params.Rs./cos(PAdots);
            Rdots = [zeros(1,length(totalLengthdots));L_now.*(1-sdots.^2).^(3/2);L_now.*(1-sdots.^2).*sdots];
            ringn = get_ring_n(Rdots);
            ringn_prob = ringn.*(totalLengthdots);
            ringnum = ringparams.ringn.*ringparams.cross.*(ringn_prob)';
            ringnum_sum = sum(ringnum,1);
            coll_length(i) = 2*sum(totalLengthdots);
            coll_num(:,i) = 2*ringnum_sum;
        end
    end
end
end