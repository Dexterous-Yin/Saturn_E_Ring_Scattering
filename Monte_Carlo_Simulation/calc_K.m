% calc K
function res = calc_K(aeq,L,params)
% aeq in deg; res in sqrt(G)*Rs
% aeq: 1xn; L: 1x1
resK = zeros(1,length(aeq));
for i=1:length(aeq)
    fun = @(x) params.Bs*1e9/L^3*sqrt(1+3*x^2)/(1-x^2)^3-params.Bs*1e9/L^3/sind(aeq(i))^2;
    sm = abs(fzero(fun,[0,0.999]));
    Bm = abs(params.Bs/L^3*sqrt(1+3*sm^2)/(1-sm^2)^3); %[T]
    deltas = sm/10000;
    sdots = -sm+deltas/2:deltas:sm-deltas/2;
    Bdots = abs(params.Bs/L^3*sqrt(1+3.*sdots.^2)./(1-sdots.^2).^3); %[T]
    slength = calc_field_length(L,sdots+deltas/2)-calc_field_length(L,sdots-deltas/2); %[Rs]
    resK(i) = sum(sqrt(Bm-Bdots).*slength);
end
res = resK*1e2;
end