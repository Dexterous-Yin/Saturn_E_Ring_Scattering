% calc mu
function res = calc_mu(W,aeq,L,params)
% W in keV, aeq in deg, res in MeV/G
W = W.*params.keV;
B = abs(params.Bs./L.^3);
mu = W.*(W+2*params.me*params.csq)./(2*params.me*params.csq.*B).*sind(aeq).^2;
res = mu./1e10./abs(params.qe);
end