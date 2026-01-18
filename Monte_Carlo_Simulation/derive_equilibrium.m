%{
This script is used to derive the equilibrium state of Saturn radiation belt electrons.
%}

% load basic settings
load('../Monte_Carlo_Simulation/Green_Function_First_Half.mat');
% bounce period for each grid
W_grid = E_grid.*params.keV;
gama_grid = W_grid./params.me/params.csq+1;
V_grid = sqrt((gama_grid.^2-1)*params.csq./gama_grid.^2);
Tb_grid = zeros(size(E_grid));
for Li = 1:length(L_grid)
    for Ki = 1:length(K_grid)
        for mui = 1:length(mu_grid)
            coll_length_temp=interp1(PA_set,squeeze(coll_length(Li,:)),PA_grid(Li,Ki,mui));
            Tb_grid(Li,Ki,mui) = coll_length_temp/V_grid(Li,Ki,mui);
        end
    end
end
disp('Tb calculation DONE.');

%% parameters
% Green function results
load('../Monte_Carlo_Simulation/Green_full.mat');
Green_full_rec = Green_full;
[K_grid_now,L_grid_now,mu_grid_now] = meshgrid(K_grid,L_grid,mu_grid);

simT = 1e4;
time_norm = 50;
n_particles = 1e6;

%%
simStep = simT/time_norm;
% Geant4 parameters
G4Emin = 120;
G4Emax = 6.2e4;
G4Enum = 300;
G4Enum_add = 29;
G4E_set = exp(interp1([1,G4Enum],[log(G4Emin),log(G4Emax)],-G4Enum_add:1:G4Enum,'linear','extrap')); % [keV]
G4E_set_edge = exp(interp1([1,G4Enum],[log(G4Emin),log(G4Emax)],-0.5-G4Enum_add:1:G4Enum+0.5,'linear','extrap')); % [keV]
G4_num = 1e6;
% modify for delted grids
for Li = 1:length(L_grid)
    for Ki = 1:length(K_grid) % Ki in
        for mui = 1:length(mu_grid) % mui in
            if del_grid(Li,Ki,mui) == 1
                Green_full(Li,Ki,mui,:,:) = 0;
                Green_full(Li,Ki,mui,Ki,mui) = n_particles;
            end
        end
    end
end
% calculate dedomega
dEdomega = zeros(size(E_grid)); % [keV*sr]
dEdsineq = zeros(size(E_grid)); % [keV*sin^2]
for Li=1:length(L_grid)
    for Ki = 1:length(K_grid)
        for mui = 1:length(mu_grid)
            domega = 2*pi*(cosd(PA_grid_edge(Li,Ki+1,mui))-cosd(PA_grid_edge(Li,Ki,mui)));
            dsineq = pi*(sind(PA_grid_edge(Li,Ki,mui))^2-sind(PA_grid_edge(Li,Ki+1,mui))^2);
            dE = (E_grid_edge(Li,Ki,mui+1)-E_grid_edge(Li,Ki,mui)+E_grid_edge(Li,Ki+1,mui+1)-E_grid_edge(Li,Ki+1,mui))/2;
            dEdomega(Li,Ki,mui) = domega*dE;
            dEdsineq(Li,Ki,mui) = dsineq*dE;
        end
    end
end
p_sq_grid = (E_grid.*params.keV).*(E_grid.*params.keV+2*params.me*params.csq)/(params.me^3*1e18*params.c^2*10/abs(params.qe)); % denominator to transfer from differential flux to psd in s^3/km^6

%% modify the Green function for phase space density
Green_psd = zeros(size(Green_full));
counts_input = zeros(size(E_grid))+n_particles;
psd_input = counts_input./dEdsineq./p_sq_grid./Tb_grid;
for Ki = 1:length(K_grid) % Ki in
    for mui = 1:length(mu_grid) % mui in
        green_now = squeeze(Green_full(:,Ki,mui,:,:));
        green_now_psd = green_now./dEdsineq./p_sq_grid./Tb_grid;
        for Li = 1:length(L_grid)
            green_ratio = squeeze(green_now_psd(Li,:,:))./psd_input(Li,Ki,mui);
            Green_psd(Li,Ki,mui,:,:) = green_ratio;
        end
    end
end

%% Dll settings
Dll_v = 3.86e-17; 
Dll_n = 13.3; 
PSD_factor = 1; 
psd_none = 0*PSD_factor;

%%
del_L = L_grid(2)-L_grid(1);
L_outbound = L_grid(end)+del_L;
E_grid_outbound = zeros(length(K_grid),length(mu_grid));
tempa = derive_PA(K_grid,L_outbound+zeros(1,length(K_grid)),mapK,mapL,mapPA);
tempaa = repmat(tempa',1,length(mu_grid));
PA_grid_outbound = tempaa;
for mui=1:length(mu_grid)
    tempE = derive_W(mu_grid(mui),squeeze(PA_grid_outbound(:,mui)),L_outbound,params);
    E_grid_outbound(:,mui) = tempE;
end

omega_b_outbound = omega_b(E_grid_outbound*params.keV,PA_grid_outbound/180*pi,L_outbound,params);
Tb_grid_outbound = 2*pi./omega_b_outbound;

sinn = 0.6;
flux_outbound = get_flux(E_grid_outbound,PA_grid_outbound,sinn);
p_sq_outbound = (E_grid_outbound.*params.keV).*(E_grid_outbound.*params.keV+2*params.me*params.csq)/(params.me^3*1e18*params.c^2*10/abs(params.qe));
psd_outbound = flux_outbound./p_sq_outbound*PSD_factor;

E_grid_edge_outbound = zeros(length(K_grid_edge),length(mu_grid_edge));
tempa = derive_PA(K_grid_edge,L_outbound+zeros(1,length(K_grid_edge)),mapK,mapL,mapPA);
tempaa = repmat(tempa',1,length(mu_grid_edge));
PA_grid_edge_outbound = tempaa;
for mui=1:length(mu_grid_edge)
    tempE = derive_W(mu_grid_edge(mui),squeeze(PA_grid_edge_outbound(:,mui)),L_outbound,params);
    E_grid_edge_outbound(:,mui) = tempE;
end

% calculate dedomega
dEdomega_outbound = zeros(size(E_grid_outbound)); % [keV*sr]
for Ki = 1:length(K_grid)
    for mui = 1:length(mu_grid)
        domega = 2*pi*(cosd(PA_grid_edge_outbound(Ki+1,mui))-cosd(PA_grid_edge_outbound(Ki,mui)));
        dE = (E_grid_edge_outbound(Ki,mui+1)-E_grid_edge_outbound(Ki,mui)+E_grid_edge_outbound(Ki+1,mui+1)-E_grid_edge_outbound(Ki+1,mui))/2;
        dEdomega_outbound(Ki,mui) = domega*dE;
    end
end
counts_outbound = flux_outbound.*dEdomega_outbound.*Tb_grid_outbound.*cosd(PA_grid_outbound);

L_inbound = L_grid(1)-del_L;
psd_inbound = psd_none;

%% construct the sparse matrix
% delete 90 deg
del_grid_now = del_grid;
del_grid_now(:,1,:) = 1;

pickid = find(del_grid_now == 0);
L_pick = L_grid_now(pickid);
K_pick = K_grid_now(pickid);
mu_pick = mu_grid_now(pickid);
E_pick = E_grid(pickid);
PA_pick = PA_grid(pickid);
Arow = length(pickid);
Acol = length(pickid);
% sparse matrix construction
sparsei = zeros(1,length(pickid));
sparsej = zeros(1,length(pickid));
sparsev = zeros(1,length(pickid));
sparseid = 1;
Bmatrix = zeros(Arow,1);
for i=1:length(pickid)
    L_now = L_pick(i);
    K_now = K_pick(i);
    mu_now = mu_pick(i);
    Lid = find(abs(L_grid-L_now)<1e-7);
    Koutid = find(abs(K_grid-K_now)<1e-7);
    muoutid = find(abs(mu_grid-mu_now)<1e-7);

    sparse_now = zeros(1,Acol);
    right_now = zeros(Arow,1);
    % first part about diffusion
    pid = find(abs(L_pick-L_now-del_L)<1e-7 & abs(K_pick-K_now)<1e-7 & abs(mu_pick-mu_now)<1e-7);
    mid = find(abs(L_pick-L_now+del_L)<1e-7 & abs(K_pick-K_now)<1e-7 & abs(mu_pick-mu_now)<1e-7);
    L_p_half = L_now+del_L/2;
    L_m_half = L_now-del_L/2;
    Dll_p_half = Dll_v*L_p_half.^Dll_n;
    Dll_m_half = Dll_v*L_m_half.^Dll_n;

    sparse_now(i) = -L_now^2/(del_L^2)*(Dll_p_half/L_p_half^2+Dll_m_half/L_m_half^2);
    if isempty(mid)
        sparse_now(pid) = L_now^2/(del_L^2)*Dll_p_half/L_p_half^2;
        right_now(i) = -L_now^2/(del_L^2)*Dll_m_half/L_m_half^2*psd_inbound;
    elseif isempty(pid)
        sparse_now(mid) = L_now^2/(del_L^2)*Dll_m_half/L_m_half^2;
        right_now(i) = -L_now^2/(del_L^2)*Dll_p_half/L_p_half^2*psd_outbound(Koutid,muoutid);
    else
        sparse_now(pid) = L_now^2/(del_L^2)*Dll_p_half/L_p_half^2;
        sparse_now(mid) = L_now^2/(del_L^2)*Dll_m_half/L_m_half^2;
    end
    sparse_now = sparse_now*simT; % Note the Green function is just PSD ratio.
    right_now = right_now*simT;

    % second part about Green function
    for m=1:length(pickid)
        if abs(L_pick(m)-L_now)<1e-7
            Kinid = find(abs(K_grid-K_pick(m))<1e-7);
            muinid = find(abs(mu_grid-mu_pick(m))<1e-7);
            sparse_now(m) = sparse_now(m)+Green_psd(Lid,Kinid,muinid,Koutid,muoutid);
        end
    end
    sparse_now(i) = sparse_now(i)-1; % delta f

    % put in full matrix
    fillid = find(sparse_now~=0);
    sparsei(sparseid:sparseid+length(fillid)-1) = zeros(1,length(fillid))+i;
    sparsej(sparseid:sparseid+length(fillid)-1) = fillid;
    sparsev(sparseid:sparseid+length(fillid)-1) = sparse_now(fillid);
    Bmatrix = Bmatrix+right_now;
    sparseid = sparseid+length(fillid);
    if mod(i,1000) == 0
        disp(i);
    end
end
% delete tail
sparsei_final = sparsei(1:sparseid-1);
sparsej_final = sparsej(1:sparseid-1);
sparsev_final = sparsev(1:sparseid-1);
Amatrix = sparse(sparsei_final,sparsej_final,sparsev_final,Arow,Acol);
% solve the system of equations
eqnres = Amatrix\Bmatrix;
%% summarize results
psd_sim = zeros(length(L_grid),length(K_grid),length(mu_grid))+psd_none;
for i=1:length(pickid)
    Lid = find(abs(L_grid-L_pick(i))<1e-7);
    Kid = find(abs(K_grid-K_pick(i))<1e-7);
    muid = find(abs(mu_grid-mu_pick(i))<1e-7);
    psd_sim(Lid,Kid,muid) = eqnres(i)/PSD_factor;
end

%% derive the distribution based on derived psd
% axis settings
del_L = L_grid(2)-L_grid(1);
L_outbound = L_grid(end)+del_L;
L_grid_full = [L_grid,L_grid(end)+del_L];
L_grid_edge_full = [L_grid_edge,L_outbound+del_L/2];
flux_dis = zeros(length(L_grid_edge_full),length(K_grid_edge));
flux_dis_line = zeros(1,length(L_grid_full));
count_dis_line = zeros(1,length(L_grid_full));
flux_dis_ratio = zeros(length(L_grid_edge_full),length(K_grid_edge));
PA_dis = zeros(length(L_grid_edge_full),length(K_grid_edge));
psd_ana = psd_sim;
flux_now = psd_ana.*p_sq_grid;
counts_now = flux_now.*dEdomega;
% combine boundary condition
psd_ana_full = zeros(length(L_grid)+1,length(K_grid),length(mu_grid));
PA_grid_full = zeros(length(L_grid)+1,length(K_grid),length(mu_grid));
PA_grid_edge_full = zeros(length(L_grid)+1,length(K_grid_edge),length(mu_grid_edge));
E_grid_edge_full = zeros(length(L_grid)+1,length(K_grid_edge),length(mu_grid_edge));
flux_full = zeros(length(L_grid)+1,length(K_grid),length(mu_grid));
counts_full = zeros(length(L_grid)+1,length(K_grid),length(mu_grid));
dEdomega_full = zeros(length(L_grid)+1,length(K_grid),length(mu_grid));

% put data in
psd_ana_full(1:length(L_grid),:,:) = psd_ana;
psd_ana_full(end,:,:) = psd_outbound/PSD_factor;
flux_full(1:length(L_grid),:,:) = flux_now;
flux_full(end,:,:) = flux_outbound;
counts_full(1:length(L_grid),:,:) = counts_now;
counts_full(end,:,:) = flux_outbound.*dEdomega_outbound;
dEdomega_full(1:length(L_grid),:,:) = dEdomega;
dEdomega_full(end,:,:) = dEdomega_outbound;
PA_grid_full(1:length(L_grid),:,:) = PA_grid;
PA_grid_full(end,:,:) = PA_grid_outbound;
PA_grid_edge_full(1:length(L_grid),:,:) = PA_grid_edge;
PA_grid_edge_full(end,:,:) = PA_grid_edge_outbound;
E_grid_edge_full(1:length(L_grid),:,:) = E_grid_edge;
E_grid_edge_full(end,:,:) = E_grid_edge_outbound;

% focused energy range
Emin = 790;
Emax = 4750; %4750;
cps_sim_full = zeros(length(L_grid)+1);
for Li = 1:length(L_grid)+1
    PA_out_edge = squeeze(PA_grid_edge_full(Li,:,1));
    PA_dis(Li,:) = PA_out_edge;
end
for Li = 1:length(L_grid)+1
    PA_out = squeeze(PA_grid_full(Li,:,1));
    PA_out_edge = squeeze(PA_grid_edge_full(Li,:,1));
    flux_out = zeros(size(PA_out_edge))+nan;
    counts_out = zeros(size(PA_out_edge))+nan;
    testEid = zeros(length(PA_grid),2);
    for PAi = 2:length(PA_out)
        counts_sum = 0;
        dEdomega_sum = 0;
        pickEid_low_before = find(squeeze(E_grid_edge_full(Li,PAi,:))<=Emin);
        if isempty(pickEid_low_before)
            pickEid_low_before=1;
        else
            pickEid_low_before=pickEid_low_before(end);
        end
        pickEid_high_before = find(squeeze(E_grid_edge_full(Li,PAi+1,:))<=Emin);
        if isempty(pickEid_high_before)
            pickEid_high_before=1;
        else
            pickEid_high_before=pickEid_high_before(end);
        end
        pickEid_before = min([pickEid_low_before,pickEid_high_before]);
        pickEid_low_after = find(squeeze(E_grid_edge_full(Li,PAi,:))>=Emax); pickEid_low_after=pickEid_low_after(1)-1;
        pickEid_high_after = find(squeeze(E_grid_edge_full(Li,PAi+1,:))>=Emax); pickEid_high_after=pickEid_high_after(1)-1;
        pickEid_after = max([pickEid_low_after,pickEid_high_after]);
        testEid(PAi,:) = [pickEid_before,pickEid_after];
        for pickEi = pickEid_before:pickEid_after
            E_edge_ul = E_grid_edge_full(Li,PAi,pickEi);
            E_edge_ur = E_grid_edge_full(Li,PAi,pickEi+1);
            E_edge_ll = E_grid_edge_full(Li,PAi+1,pickEi);
            E_edge_lr = E_grid_edge_full(Li,PAi+1,pickEi+1);
            ratio = 0;
            if E_edge_ul>=Emin && E_edge_ll>=Emin && E_edge_ur<=Emax && E_edge_lr<=Emax
                ratio = 1;
            elseif E_edge_ul<=Emin && E_edge_ll<=Emin && E_edge_ur>=Emin && E_edge_lr>=Emin
                ratio = (E_edge_ur-Emin+E_edge_lr-Emin)/(E_edge_ur-E_edge_ul+E_edge_lr-E_edge_ll);
            elseif E_edge_ul<Emin && E_edge_ll<Emin && E_edge_ur<Emin && E_edge_lr>Emin
                ratio = (E_edge_lr-Emin)*(E_edge_lr-Emin)/(E_edge_lr-E_edge_ur)/(E_edge_ur-E_edge_ul+E_edge_lr-E_edge_ll);
            elseif E_edge_ul<Emin && E_edge_ll>Emin && E_edge_ur<=Emax && E_edge_lr<=Emax
                ratio = 1-(Emin-E_edge_ul)*(Emin-E_edge_ul)/(E_edge_ll-E_edge_ul)/(E_edge_ur-E_edge_ul+E_edge_lr-E_edge_ll);
            elseif E_edge_ul<=Emax && E_edge_ll<=Emax && E_edge_ur>=Emax && E_edge_lr>=Emax
                ratio = (Emax-E_edge_ul+Emax-E_edge_ll)/(E_edge_ur-E_edge_ul+E_edge_lr-E_edge_ll);
            elseif E_edge_ul<Emax && E_edge_ll>Emax && E_edge_ur>Emax && E_edge_lr>Emax
                ratio = (Emax-E_edge_ul)*(Emax-E_edge_ul)/(E_edge_ll-E_edge_ul)/(E_edge_ur-E_edge_ul+E_edge_lr-E_edge_ll);
            elseif E_edge_ul<Emax && E_edge_ll<Emax && E_edge_ur<Emax && E_edge_lr>Emax
                ratio = 1-(E_edge_lr-Emax)*(E_edge_lr-Emax)/(E_edge_lr-E_edge_ur)/(E_edge_ur-E_edge_ul+E_edge_lr-E_edge_ll);
            elseif E_edge_ul<Emax && E_edge_ll>Emax && E_edge_ur<Emax && E_edge_lr>Emax
                ratio = (Emax-E_edge_ul)*(Emax-E_edge_ul)/(E_edge_ll-E_edge_ul)/(E_edge_ur-E_edge_ul+E_edge_lr-E_edge_ll)-(Emax-E_edge_ur)*(Emax-E_edge_ur)/(E_edge_lr-E_edge_ur)/(E_edge_ur-E_edge_ul+E_edge_lr-E_edge_ll);
            else
                disp('Note Here.');
            end
            counts_sum = counts_sum+counts_full(Li,PAi,pickEi)*ratio;
            dEdomega_sum = dEdomega_sum+dEdomega_full(Li,PAi,pickEi)*ratio;
        end
        dEdomega_total = 2*pi*(cosd(PA_grid_edge_full(Li,PAi+1,1))-cosd(PA_grid_edge_full(Li,PAi,1)))*(Emax-Emin);
        flux_out(PAi) = counts_sum/dEdomega_sum;
        counts_out(PAi) = counts_sum;
    end
    % flux_ave = sum(flux_out(1:length(PA_out)).*(cosd(PA_out_edge(2:end))-cosd(PA_out_edge(1:end-1))),'omitnan');
    temp_flux = flux_out(1:length(PA_out));
    temp_counts = counts_out(1:length(PA_out));
    % temp_del_PA = (PA_out_edge(1:end-1))-(PA_out_edge(2:end));
    temp_del_PA = cosd(PA_out_edge(2:end))-cosd(PA_out_edge(1:end-1));
    temp_del_PA(isnan(temp_flux)) = nan;
    flux_ave = sum(temp_flux.*temp_del_PA,'omitnan')/sum(temp_del_PA,'omitnan');
    count_ave = sum(temp_counts.*temp_del_PA,'omitnan')/sum(temp_del_PA,'omitnan');
    % flux_ave = mean(flux_out(1:length(PA_out)));
    flux_out_ratio = flux_out./flux_ave;
    flux_dis(Li,:) = flux_out;
    flux_dis_ratio(Li,:) = flux_out_ratio;
    flux_dis_line(Li) = flux_ave;
    count_dis_line(Li) = count_ave;
end

%% save
save('./Data/Equilibrium_result.mat','L_grid_edge_full','K_grid_edge','PA_dis','flux_dis_ratio','L_grid_full','flux_dis_line');

%% functions

function res = get_flux(E,PA,n)
f0 = 138.60386;
E0 = 100;
gamma = -1.036261;
Ec = 2851.9122;
KT = 520.18468;
res = f0.*(E./E0).^gamma./(1+exp((E-Ec)./KT)).*sind(PA).^n;
end