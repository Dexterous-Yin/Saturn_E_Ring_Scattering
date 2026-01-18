[params, ~, ~] = set_parameters('e');
L = 2.3:0.01:5.1;
PA = 2:0.01:90;
[mapPA,mapL] = meshgrid(PA,L);
mapK = zeros(length(L),length(PA));
tempnum = length(PA);
parfor Li = 1:length(L)
    for PAi = 1:tempnum
        tempK = calc_K(PA(PAi),L(Li),params);
        mapK(Li,PAi) = tempK;
    end
    disp(Li);
end

save('./Kmap.mat','mapK','mapL','mapPA');


% calc K
function res = calc_K(aeq,L,params)
% aeq in deg % res in sqrt(G)*Rs
fun = @(x) params.Bs*1e9/L^3*sqrt(1+3*x^2)/(1-x^2)^3-params.Bs*1e9/L^3/sind(aeq)^2;
sm = abs(fzero(fun,[0,0.999]));
Bm = abs(params.Bs/L^3*sqrt(1+3*sm^2)/(1-sm^2)^3); %[T]
deltas = sm/10000;
sdots = -sm+deltas/2:deltas:sm-deltas/2;
Bdots = abs(params.Bs/L^3*sqrt(1+3.*sdots.^2)./(1-sdots.^2).^3); %[T]
slength = calc_field_length(L,sdots+deltas/2)-calc_field_length(L,sdots-deltas/2); %[Rs]
K = sum(sqrt(Bm-Bdots).*slength);
res = K*1e2;
end
% field line length
function res = calc_field_length(L,s)
res = L.*(asinh(sqrt(3).*s)./(2*sqrt(3))+s.*sqrt(1+3*s.^2)./2);
end

function [params, scaling, Charge]= set_parameters(species)
    % fundamental parameters
    params.qe = -1.6021766e-19;
    params.keV = abs(1e3*params.qe);
    params.me = 9.1093837e-31;
    params.c = 299792458;
    params.csq = params.c^2;

    % the magnetic field
    params.Bs = -21e-6; % equatorial magnetic field at L=1
    params.Rs = 60268e3; %[m]
    
    % particle species
    if nargin == 0
        Mass = 1; Charge = -1;
    elseif species == 'e'
        Mass = 1; Charge = -1;
    end
    
    % calculates scaling coefficients for the variables    
    scaling.B = Mass*params.me*params.c/abs(Charge)/abs(params.qe)/params.Rs/params.Bs; % [UnitB/Bs]
    scaling.E = Mass*params.me*params.csq/abs(Charge)/abs(params.qe)/params.Rs/(1e-3); % [UnitE/(mV/m)]
    scaling.r = 1.0; % [UnitR/Rs]
    scaling.Energy = Mass*params.me*params.csq/params.keV; % [UnitE/keV]
    scaling.Time = params.Rs/params.c; % [UnitT]
    scaling.Angle = 180/pi; % [UnitAngle]
    scaling.p_par = 1; % [UnitP/(m0*c)]
    
end