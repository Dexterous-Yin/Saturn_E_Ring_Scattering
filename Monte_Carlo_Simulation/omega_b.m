% bounce omega
function res = omega_b(W,aeq,L,params)
gamma = 1+W./(params.me*params.c^2);
v = params.c*sqrt(1-1./gamma.^2);
% res = pi.*v./(2.*L.*params.Rs.*(1.3-0.56.*sin(aeq)));
res = pi.*v./(2.*L.*params.Rs.*(1.3802-0.3198.*(sin(aeq)+sin(aeq).^(1/2))));
end