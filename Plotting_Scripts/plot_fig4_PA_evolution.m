%{
This script is used to simulate the temporal evolution of pitch-angle distribution.
%}
% load matplotlib.pyplot
pyplt = py.importlib.import_module('matplotlib.pyplot');
cmap = pyplt.get_cmap('viridis');
cmap_vals = cmap(linspace(0, 1, 256));
cmap_mat = double(py.array.array('d', py.numpy.nditer(cmap_vals)));
cmap_mat = reshape(cmap_mat, [4, 256])';
colormapData = cmap_mat(:,1:3);
colormap_viridis = colormapData; %turbo;

%% load basic settings
load('../Monte_Carlo_Simulation/Data/Green_Function_First_Half.mat');
% bounce period for each grid
[K_grid_now,L_grid_now,mu_grid_now] = meshgrid(K_grid,L_grid,mu_grid);
W_grid = E_grid.*params.keV;
gama_grid = W_grid./params.me/params.csq+1;
V_grid = sqrt((gama_grid.^2-1)*params.csq./gama_grid.^2);
Tb_grid = zeros(size(E_grid));
for Li = 17
    for Ki = 1:length(K_grid)
        for mui = 1:length(mu_grid)
            coll_length_temp=interp1(PA_set,squeeze(coll_length(Li,:)),PA_grid(Li,Ki,mui));
            Tb_grid(Li,Ki,mui) = coll_length_temp/V_grid(Li,Ki,mui);
        end
    end
end
disp('Tb calculation DONE.');

%% parameters
n_particles = 1e6;
time_norm = 50; %[s]
simT = 1e4;
% Green function results
load(['../Monte_Carlo_Simulation/Data/Green_full.mat']);

simStep = simT/time_norm;
% modify for delted grids
for Li = 17 %1:length(L_grid)
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
for Li = 17 %1:length(L_grid)
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
        for Li = 17 %1:length(L_grid)
            green_ratio = squeeze(green_now_psd(Li,:,:))./psd_input(Li,Ki,mui);
            Green_psd(Li,Ki,mui,:,:) = green_ratio;
        end
    end
end

%% lifetime
Green_lifetime = zeros(length(L_grid),length(K_grid),length(mu_grid));
Green_num = zeros(length(L_grid),length(K_grid),length(mu_grid));
for Li = 17 
    for Ki = 1:length(K_grid) % Ki in
        for mui = 1:length(mu_grid) % mui in
            green_ratio = Green_full(Li,Ki,mui,Ki,mui)/n_particles;
            Green_num(Li,Ki,mui) = Green_full(Li,Ki,mui,Ki,mui);
            Green_lifetime(Li,Ki,mui) = -simT/log(green_ratio);
            if green_ratio == 1
                Green_lifetime(Li,Ki,mui) = nan;
            end
        end
    end
end

%%
Lpick = 17;
flux_ini = zeros(length(K_grid),length(mu_grid));
counts_ini = zeros(length(K_grid),length(mu_grid));
psd_ini = zeros(length(K_grid),length(mu_grid));

sinn = 0.6;
PSD_factor = 1;
for Ki = 1:length(K_grid) % Ki in
    for mui = 1:length(mu_grid) % mui in
        flux_ini(Ki,mui) = get_flux(E_grid(Lpick,Ki,mui),PA_grid(Lpick,Ki,mui),sinn);
        counts_ini(Ki,mui) = flux_ini(Ki,mui)*dEdomega(Lpick,Ki,mui);
        psd_ini(Ki,mui) = flux_ini(Ki,mui)/p_sq_grid(Lpick,Ki,mui);
    end
end

del_grid_now = del_grid;
del_grid_now(:,1,:) = 1;
del_grid_pick = squeeze(del_grid_now(Lpick,:,:));

flux_ini(del_grid_pick==1) = 0;
counts_ini(del_grid_pick==1) = 0;
psd_ini(del_grid_pick==1) = 0;
flux_ini = flux_ini*PSD_factor;
counts_ini = counts_ini*PSD_factor;
psd_ini = psd_ini*PSD_factor;

Green_psd_pick = squeeze(Green_psd(Lpick,:,:,:,:));

%%
disp('Simulation for temporal evolution begins...');
steps = 0:1:5000;
counts = zeros(length(steps),length(K_grid),length(mu_grid));
counts(1,:,:) = counts_ini;
psd_sim = zeros(length(steps),length(K_grid),length(mu_grid));
psd_sim(1,:,:) = psd_ini;
flux_sim = zeros(length(steps),length(K_grid),length(mu_grid));
flux_sim(1,:,:) = flux_ini;
for si=1:length(steps)-1
    tempout_psd = zeros(length(K_grid),length(mu_grid));
    f_l = squeeze(psd_sim(si,:,:));
    for Ki = 1:length(K_grid) % Ki in
        for mui = 1:length(mu_grid) % mui in
            tempout_psd(Ki,mui) = sum(f_l.*squeeze(Green_psd_pick(:,:,Ki,mui)),'all','omitnan');
        end
    end
    tempout_psd(del_grid_pick == 1) = 0;
    tempout_psd(~isfinite(tempout_psd)) = nan;
    % counts(si+1,:,:,:) = tempout;
    psd_sim(si+1,:,:) = tempout_psd;
    flux_sim(si+1,:,:) = tempout_psd.*squeeze(p_sq_grid(Lpick,:,:));
    counts(si+1,:,:) = tempout_psd.*squeeze(p_sq_grid(Lpick,:,:)).*squeeze(dEdomega(Lpick,:,:));
    if mod(si,10)==0
        disp([num2str(si),'/',num2str(5000)]);
    end
end
psd_sim = psd_sim/PSD_factor;
flux_sim = flux_sim/PSD_factor;
counts = counts/PSD_factor;

%% original data
Emin = 790;
Emax = 4750;

PA_grid_now = squeeze(PA_grid(Lpick,:,:));
E_grid_now = squeeze(E_grid(Lpick,:,:));
PA_out = PA_grid_now(:,1);
PA_out_edge = squeeze(PA_grid_edge(Lpick,:,1));
time_out = steps.*simT;
time_out_edge = [-0.5:1:max(steps)+1].*simT;
flux_out = zeros(length(PA_out_edge),length(time_out_edge))+nan;
flux_out_ratio = zeros(length(PA_out_edge),length(time_out_edge))+nan;
flux_out_line = zeros(1,length(time_out_edge))+nan;
rec = [];
dedomegarec = [];

peak_pa_plot = zeros(1,length(time_out));

for si=1:length(steps)
    for PAi=2:length(PA_out)
        counts_sum = 0;
        dEdomega_sum = 0;
        pickPAid = find(PA_grid_now(:,1)>=PA_out_edge(PAi+1) & PA_grid_now(:,1)<=PA_out_edge(PAi));
        for pickPAi = pickPAid(1):pickPAid(end)
            pickEid_low_before = find(squeeze(E_grid_edge(Lpick,pickPAi,:))<=Emin);
            if isempty(pickEid_low_before)
                pickEid_low_before=1;
            else
                pickEid_low_before=pickEid_low_before(end);
            end
            pickEid_high_before = find(squeeze(E_grid_edge(Lpick,pickPAi+1,:))<=Emin);
            if isempty(pickEid_high_before)
                pickEid_high_before=1;
            else
                pickEid_high_before=pickEid_high_before(end);
            end
            pickEid_before = min([pickEid_low_before,pickEid_high_before]);
            pickEid_low_after = find(squeeze(E_grid_edge(Lpick,pickPAi,:))>=Emax); pickEid_low_after=pickEid_low_after(1)-1;
            pickEid_high_after = find(squeeze(E_grid_edge(Lpick,pickPAi+1,:))>=Emax); pickEid_high_after=pickEid_high_after(1)-1;
            pickEid_after = max([pickEid_low_after,pickEid_high_after]);
            for pickEi = pickEid_before:pickEid_after
                E_edge_ul = E_grid_edge(Lpick,pickPAi,pickEi);
                E_edge_ur = E_grid_edge(Lpick,pickPAi,pickEi+1);
                E_edge_ll = E_grid_edge(Lpick,pickPAi+1,pickEi);
                E_edge_lr = E_grid_edge(Lpick,pickPAi+1,pickEi+1);
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
                counts_sum = counts_sum+counts(si,pickPAi,pickEi)*ratio;
                dEdomega_sum = dEdomega_sum+dEdomega(Lpick,pickPAi,pickEi)*ratio;
            end
        end
        dEdomega_total = 2*pi*(cosd(PA_grid_edge(Lpick,PAi+1,1))-cosd(PA_grid_edge(Lpick,PAi,1)))*(Emax-Emin);
        flux_out(PAi,si) = counts_sum/dEdomega_sum;
    end
    temp_flux = flux_out(1:length(PA_out),si)';
    temp_del_PA = cosd(PA_out_edge(2:end))-cosd(PA_out_edge(1:end-1));
    temp_del_PA(isnan(temp_flux)) = nan;
    flux_ave = sum(temp_flux.*temp_del_PA,'omitnan')/sum(temp_del_PA,'omitnan');
    flux_out_ratio(:,si) = flux_out(:,si)./flux_ave;

    if mod(si,1000) == 0
        disp(si);
    end
end

%% combine with life time

fig = figure('Color',[1 1 1]);
layout = tiledlayout(2,1, 'TileSpacing', 'tight', 'Padding', 'compact');

ax1 = nexttile;
flux_out_ratio(flux_out_ratio==0)=nan;
flux_out_ratio_plot = flux_out_ratio(1:length(PA_out),1:length(time_out));

pco_nan = zeros(size(flux_out_ratio_plot));
pco_nan(1,:) = nan;
flux_out_ratio_plot(1,:) = nan;
PA_out_plot = PA_out;
pco = imagesc(time_out,PA_out_plot(1:end),flux_out_ratio_plot(1:end,:));
set(pco, 'AlphaData', ~isnan(pco_nan));

% pco = pcolor(time_out_edge,PA_out_edge,flux_out_ratio);
% set(pco,'EdgeColor','none');
% hold on;
colormap(gca,colormap_viridis);
c = colorbar('FontSize',14);
c.LineWidth = 1.5;
c.Label.String = 'Flux/omniFlux';
c.Label.FontSize = 16;
c.TickLength = 0.015;
% clim([0,2]);
clim([0.2,1.23]);
ylim([10,90]);
xlim([0,5e7]);
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize', 14, 'Layer','top','Color',[0.8,0.8,0.8]);
set(gca,'XTick',0:1e7:5e7);
xlabel('Time (s)','FontSize',16);
ylabel('PA (deg)','FontSize',16);
set(gca,'TickDir','out','Box','on','YDir','normal');
% title(['[',num2str(Emin),',',num2str(Emax),' keV]']);

% life time
caseL = 17;
tempGreen = squeeze(Green_lifetime(caseL,:,:));

pick_PA_grid_edge = squeeze(PA_grid_edge(caseL,:,:));
pick_E_grid_edge = squeeze(E_grid_edge(caseL,:,:));
pick_Green_grid = zeros(length(K_grid_edge),length(mu_grid_edge))+nan;
pick_Green_grid(1:length(K_grid),1:length(mu_grid)) = tempGreen;
pick_Green_grid(1,:) = nan;

pick_E_grid = squeeze(E_grid(caseL,:,:));
pick_PA_grid = squeeze(PA_grid(caseL,:,:));
pick_mu_grid = squeeze(mu_grid_now(caseL,:,:));
pick_E_range = [400,10e3];
pick_PA_range = [20,88];
pick_range_id = find(pick_E_grid>=pick_E_range(1) & pick_E_grid<=pick_E_range(2) & pick_PA_grid>=pick_PA_range(1) & pick_PA_grid<pick_PA_range(2));
pick_grid_in = tempGreen(pick_range_id);
pick_lifetime_range = [min(pick_grid_in,[],'all'),max(pick_grid_in,[],'all')];
pick_lifetime_range(2) = min([pick_lifetime_range(2),5e8]);
pick_lifetime_range(1) = min([pick_lifetime_range(1),1e7]);

pick_mu_grid_in = pick_mu_grid(pick_range_id);
pick_mu_range = [min(pick_mu_grid_in,[],'all'),max(pick_mu_grid_in,[],'all')];

pick_E_range = [400,10e3];
pick_PA_range = [10,90];

nexttile;
pco2 = pcolor(pick_E_grid_edge,pick_PA_grid_edge,pick_Green_grid);
hold on;
set(pco2,'EdgeColor','none','FaceColor','flat');
set(gca,'XScale','log','ColorScale','log');
colormap(gca,colormap_viridis(end:-1:1,:));
clim(pick_lifetime_range);
c = colorbar('FontSize',14);
c.Label.FontSize = 16;
c.LineWidth = 1.5;
c.TickLength = 0.015;
% c.Label.String = 'LifeTime (s)';
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize', 14,'Color',[0.8,0.8,0.8], 'Layer','top');
set(gca,'TickDir','out');
xlabel('E (keV)','FontSize',16);
ylabel('PA (deg)','FontSize',16);
xlim(pick_E_range);
% sgtitle('Green Function','FontSize',16);
ylim(pick_PA_range);

[conm,conc] = contour(pick_E_grid_edge,pick_PA_grid_edge,pick_Green_grid,[1e5,3e5,1e6,3e6],'w--','LineWidth',2);


% window control
% screen size
set(0, 'Units', 'pixels');
screenSize = get(0, 'ScreenSize');
screenHeight = screenSize(4);

% set the height of figure to the screen height
fig.Units = 'pixels';
fig.Position(4) = 710;
fig.Position(3) = 480; %670;
fig.Position(1) = (screenSize(3) - fig.Position(3)) / 2;
fig.Position(2) = (screenSize(4) - fig.Position(4)) / 2;

%%

function res = get_flux(E,PA,n)
f0 = 138.60386;
E0 = 100;
gamma = -1.036261;
Ec = 2851.9122;
KT = 520.18468;
res = f0.*(E./E0).^gamma./(1+exp((E-Ec)./KT)).*sind(PA).^n;
end
