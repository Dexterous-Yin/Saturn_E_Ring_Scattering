% derive W from mu, L, and PA
function W = derive_W(mu,PA,L,params)
% mu in MeV/G, PA in deg, W in keV
Blocal = abs(params.Bs./L^3); %[T]
psq = mu*1e10*abs(params.qe)*2*params.me*Blocal./sind(PA).^2;
gama = sqrt(psq./(params.me*params.c)^2+1);
W = (gama-1)*params.me*params.csq/params.keV;
end