%{
This script is used to calculate the Green function at each input grid (L,K,mu).
Modify "for Li = 34" in Line 39 to loop through all L-shells.
%}
% parameters
load('./Green_Function_First_Half.mat');
% modify the following path to the folder containing Geant4 files for each energy and grain size
datapath = '/to_be_modified';

% rand seed
rng('shuffle', 'threefry');

%% calculate Green function
% Green function
Green = zeros(length(L_grid),length(K_grid),length(mu_grid),length(K_grid),length(mu_grid)); % L, Kin, muin, Kout, muout
Collnum = zeros(length(L_grid),length(K_grid),length(mu_grid));
Daa = zeros(length(L_grid),length(K_grid),length(mu_grid));
mean_delPA = zeros(length(L_grid),length(K_grid),length(mu_grid));
Corr_rec = zeros(length(L_grid),length(K_grid),length(mu_grid));
% parameters
n_particles = 1e6;
time_norm = 50; %[s]
simT = 1e4;
simStep = simT/time_norm;
% Geant4 parameters
G4Emin = 120;
G4Emax = 6.2e4;
G4Enum = 300;
G4E_set = exp(interp1([1,G4Enum],[log(G4Emin),log(G4Emax)],1:1:G4Enum)); % [keV]
G4E_set_edge = exp(interp1([1,G4Enum],[log(G4Emin),log(G4Emax)],0.5:1:G4Enum+0.5,'linear','extrap')); % [keV]
G4_num = 1e6;
% parfor
mypar = parpool(64);
% construct for parfor
N_K = length(K_grid);
N_mu = length(mu_grid);
N_states = N_K * N_mu;
% simulation
for Li = 34
    Green_now = zeros(length(K_grid),length(mu_grid),length(K_grid),length(mu_grid)); % Kin, muin, Kout, muout
    tic
    for Ki = 1:length(K_grid) % Ki in
        for mui = 1:length(mu_grid) % mui in
            disp(['---------------------------(',num2str(Li),' , ',num2str(Ki),' , ',num2str(mui),')---------------------------']);
            if del_grid(Li,Ki,mui)==1
                continue
            end
            E_ini = E_grid(Li,Ki,mui);
            PA_ini = PA_grid(Li,Ki,mui);
            E_particles = zeros(1,n_particles)+E_ini;
            PA_particles = zeros(1,n_particles)+PA_ini;
            Coll_num_particles = zeros(1,n_particles);
            % setup recorder
            recordT = 100;
            recordStep = simT/recordT+1;
            record_E_particles = zeros(recordStep,n_particles)+nan;
            record_PA_particles = zeros(recordStep,n_particles)+nan;
            % initial recorder
            recordi = 1;
            record_E_particles(recordi,:) = E_particles;
            record_PA_particles(recordi,:) = PA_particles;
            for stepi=1:simStep
                E_bound_id = discretize(E_particles,G4E_set_edge);
                E_particles_next = zeros(1,n_particles)+nan;
                PA_particles_next = zeros(1,n_particles)+nan;
                coll_length_temp = interp1(PA_set,coll_length(Li,:),PA_particles);
                W_temp = E_particles.*params.keV;
                gama_temp = W_temp./params.me/params.csq+1;
                v_temp = sqrt((gama_temp.^2-1)*params.csq./gama_temp.^2);
                coll_time_temp = coll_length_temp./v_temp;
                coll_num_temp = interp1(PA_set,squeeze(coll_num(Li,:,:)),PA_particles);
                coll_num_temp = coll_num_temp.*(time_norm./coll_time_temp)';
                parfor particlei=1:n_particles
                    if isnan(coll_length_temp(particlei)) || isnan(E_bound_id(particlei))
                        E_temp = nan; PA_temp = nan;
                    else
                        coll_num_now = coll_num_temp(particlei,:);
                        % Check if any value in coll_num_now exceeds 1
                        if any(coll_num_now > 1)
                            error('Error: coll_num_now exceeds 1. Maximum value: %.4f', max(coll_num_now));
                        end
                        rand_num = rand(1,length(ringparams.size));
                        pickid = find(rand_num<=coll_num_now);
                        E_temp = E_particles(particlei);
                        PA_temp = PA_particles(particlei); % equatorial pitch angle
                        for ri=1:length(pickid)
                            Ecollid = discretize(E_temp,G4E_set_edge);
                            if isnan(Ecollid)
                                E_temp = nan; PA_temp = nan;
                            else
                                sizeid = pickid(ri);
                                G4_res = load([datapath,num2str(G4E_set(Ecollid),'%.4f'),'_',num2str(ringparams.size(sizeid)*2,'%.2f'),'_result.mat']);
                                changeid = randi(G4_num);
                                if changeid >length(G4_res.E)
                                    E_temp = nan; PA_temp = nan;
                                else
                                    % pitch angle variation
                                    phi = 2*pi*rand(); % the random azimuthal angle change around the initial momentum direction
                                    changeangle = G4_res.PA(changeid);
                                    coll_lat = calc_coll_position(L_grid(Li),PA_temp,params);
                                    Bnow = sqrt(1+3*sind(coll_lat).^2)./cosd(coll_lat).^6;
                                    PAnow_sin = sqrt(Bnow*sind(PA_temp)^2);
                                    if PAnow_sin>1
                                        PAnow_sin = 1;
                                    end
                                    PAnow = asin(PAnow_sin); % local pitch angle
                                    PAout = calc_newPA(changeangle/180*pi,PAnow,phi);
                                    Coll_num_particles(particlei) = Coll_num_particles(particlei)+1;
                                    PA_temp = asin(sqrt(sin(PAout)^2/Bnow))/pi*180; % back to equatorial pitch angle
                                    E_temp = E_temp+G4_res.E(changeid)-G4E_set(Ecollid);
                                end
                            end
                        end
                    end
                    E_particles_next(particlei) = E_temp;
                    PA_particles_next(particlei) = PA_temp;
                end
                E_particles = E_particles_next;
                PA_particles = PA_particles_next;
                if mod(stepi,recordT/time_norm)==0
                    % disp(stepi);
                    recordi = recordi+1;
                    record_E_particles(recordi,:) = E_particles;
                    record_PA_particles(recordi,:) = PA_particles;
                end
            end
            mu_particles = calc_mu(record_E_particles(end,:),record_PA_particles(end,:),L_grid(Li),params);
            K_particles = derive_K(record_PA_particles(end,:),L_grid(Li)+zeros(1,n_particles),mapK,mapL,mapPA); %calc_K(PA_particles,L_grid(Li),params);
            mu_particles_id = discretize(mu_particles,mu_grid);
            mu_mod_id = find(isnan(mu_particles_id) & abs(mu_particles-mu_grid(1))<1e-7);
            mu_particles_id(mu_mod_id) = 1;
            mu_particles(mu_mod_id) = mu_grid(1);
            K_particles_id = discretize(K_particles,K_grid);
            % CIC process
            for particlei= 1:n_particles
                if isfinite(mu_particles_id(particlei)) && isfinite(K_particles_id(particlei))
                    % CIC based on log(mu) and PA
                    cic_mu_now = mu_particles(particlei);
                    cic_mu_left = mu_grid(mu_particles_id(particlei));
                    cic_mu_right = mu_grid(mu_particles_id(particlei)+1);
                    cic_PA_now = record_PA_particles(end,particlei);
                    cic_PA_low = PA_grid(Li,K_particles_id(particlei),1);
                    cic_PA_up = PA_grid(Li,K_particles_id(particlei)+1,1);
                    dx = (log(cic_mu_now)-log(cic_mu_left))/(log(cic_mu_right)-log(cic_mu_left));
                    dy = (cic_PA_low-cic_PA_now)/(cic_PA_low-cic_PA_up);
                    % values
                    Green_now(Ki,mui,K_particles_id(particlei),mu_particles_id(particlei)) = Green_now(Ki,mui,K_particles_id(particlei),mu_particles_id(particlei))+(1-dx)*(1-dy);
                    Green_now(Ki,mui,K_particles_id(particlei),mu_particles_id(particlei)+1) = Green_now(Ki,mui,K_particles_id(particlei),mu_particles_id(particlei)+1)+dx*(1-dy);
                    Green_now(Ki,mui,K_particles_id(particlei)+1,mu_particles_id(particlei)) = Green_now(Ki,mui,K_particles_id(particlei)+1,mu_particles_id(particlei))+(1-dx)*dy;
                    Green_now(Ki,mui,K_particles_id(particlei)+1,mu_particles_id(particlei)+1) = Green_now(Ki,mui,K_particles_id(particlei)+1,mu_particles_id(particlei)+1)+dx*dy;
                end
            end
            Collnum(Li,Ki,mui) = mean(Coll_num_particles);    
        end
    end
    duration = toc;
    Green(Li,:,:,:,:) = Green_now;
    save(['./Green_Function_L_',num2str(L_grid(Li),'%.2f'),'_CIC_',datestr(now,'mmddHHMMSS'),'.mat'],'duration','Green','n_particles','simT','time_norm','recordT','Collnum','G4_num','-v7.3');
end
delete(mypar);