% calc collision occurance position
function coll_lat = calc_coll_position(L,PA,params)
coll_lat = zeros(1,length(PA))+nan;
L_now = L;
lc_eq = asin(sqrt(1/sqrt(4*L_now^6-3*L_now^5)));
pa_now = PA;
if pa_now>=lc_eq/pi*180
    if abs(pa_now-90)<=1e-3
        coll_lat = 0;
    else
        F = @(x) params.Bs*1e9/L_now^3*sqrt(1+3*x^2)/(1-x^2)^3-params.Bs*1e9/L_now^3/sind(pa_now)^2;
        sm = abs(fzero(F,[0,0.999]));
        deltas = sm/2000;
        sdots = -sm+deltas/2:deltas:sm-deltas/2;
        Bdots = sqrt(1+3.*sdots.^2)./(1-sdots.^2).^3;
        PAdots = asin(sqrt(Bdots.*sind(pa_now)^2));
        Lengthdots = calc_field_length(L_now,sdots+deltas/2)-calc_field_length(L_now,sdots-deltas/2);
        totalLengthdots = Lengthdots*params.Rs./cos(PAdots);
        Rdots = [zeros(1,length(totalLengthdots));L_now.*(1-sdots.^2).^(3/2);L_now.*(1-sdots.^2).*sdots];
        ringn = get_ring_n(Rdots);
        ringn_prob = ringn.*(totalLengthdots);
        ringn_prob_sum = cumsum(ringn_prob);
        ringn_prob_sum = ringn_prob_sum./ringn_prob_sum(end);
        ringn_prob_sum = [0,ringn_prob_sum];
        lat_dots = asin([-sm,sdots+deltas/2])/pi*180;
        temprand = rand();
        [ring_prob_sum_unique,unique_id,~] = unique(ringn_prob_sum);
        coll_lat = interp1(ring_prob_sum_unique,lat_dots(unique_id),temprand);
    end
end
end