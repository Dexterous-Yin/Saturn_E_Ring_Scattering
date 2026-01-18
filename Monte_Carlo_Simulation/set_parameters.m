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
    params.Omegas = 2*pi/10.7/60/60;
    
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



