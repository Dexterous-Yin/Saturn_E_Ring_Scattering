%{
This script is used to plot the equilibrium state derived in ../Monte_Carlo_Simulation/derive_equilibrium.m
%}
% load matplotlib.pyplot
pyplt = py.importlib.import_module('matplotlib.pyplot');
cmap = pyplt.get_cmap('viridis');
cmap_vals = cmap(linspace(0, 1, 256));
cmap_mat = double(py.array.array('d', py.numpy.nditer(cmap_vals)));
cmap_mat = reshape(cmap_mat, [4, 256])';
colormapData = cmap_mat(:,1:3);
colormap_viridis = colormapData;

% load basic settings
load('../Monte_Carlo_Simulation/Data/Green_Function_First_Half.mat');
load('../Monte_Carlo_Simulation/Data/Equilibrium_result.mat');

%% compare with observations
% load observation data
obs_filename = '../Monte_Carlo_Simulation/Data/statistics_cassini.mat';
load(obs_filename);
obs_L_grid = double(obs_data.l_grid);
obs_L_grid_bound = double(obs_data.l_grid_bound);
obs_PA_grid = double(obs_data.pa_grid);
obs_PA_grid_bound = double(obs_data.pa_grid_bound);
obs_PA_grid_bound_rad = double(obs_PA_grid_bound)/180*pi;
obs_PA_grid_delcos = cosd(obs_PA_grid_bound(1:end-1))-cosd(obs_PA_grid_bound(2:end));

% HET CPS
hete_E_low = 790/1e3;
hete_E_high = 4750/1e3;
hete_cps_now = double(obs_data.hete_cps_median);
hete_cps_ratio_now = double(obs_data.hete_cps_ratio_median);
hete_cps_num_now = double(obs_data.hete_cps_num);
hete_cps_sem_now = double(obs_data.hete_cps_sem);
% delete
delid = find(hete_cps_num_now<3);
hete_cps_ratio_now(delid) = nan;
hete_cps_sem_now(delid) = nan;

hete_cps_now_line = zeros(1,length(obs_L_grid));
hete_cps_sem_line = zeros(1,length(obs_L_grid));
hete_cps_calcratio = zeros(size(hete_cps_ratio_now));
for Li = 1:length(obs_L_grid)
    temp_dis = hete_cps_now(Li,:);
    temp_ave = sum(temp_dis.*obs_PA_grid_delcos,'omitnan')/sum(obs_PA_grid_delcos(isfinite(temp_dis)));
    hete_cps_calcratio(Li,:) = temp_dis./temp_ave;
    hete_cps_now_line(Li) = temp_ave;
end
for Li = 1:length(obs_L_grid)
    temp_sem = hete_cps_sem_now(Li,:);
    temp_w = obs_PA_grid_delcos;
    temp_ave = sqrt(sum((temp_sem.^2).*(temp_w.^2),'omitnan'))/sum(temp_w(isfinite(temp_sem)));
    hete_cps_sem_line(Li) = temp_ave;
end
% plot variables
hete_cps_plot = zeros(length(obs_PA_grid_bound),length(obs_L_grid_bound))+nan;
hete_cps_ratio_plot = zeros(length(obs_PA_grid_bound),length(obs_L_grid_bound))+nan;
hete_cps_calcratio_plot = zeros(length(obs_PA_grid_bound),length(obs_L_grid_bound))+nan;
hete_cps_plot(1:length(obs_PA_grid),1:length(obs_L_grid)) = hete_cps_now';
hete_cps_ratio_plot(1:length(obs_PA_grid),1:length(obs_L_grid)) = hete_cps_ratio_now';
hete_cps_calcratio_plot(1:length(obs_PA_grid),1:length(obs_L_grid)) = hete_cps_calcratio';

% plot
fig = figure('Color',[1 1 1]);
layout = tiledlayout(3,1, 'TileSpacing', 'tight', 'Padding', 'compact');
xrange = [2.4,5.1];
yrange = [15,90];

clim_range = [0.2,1.3];
nexttile
hete_ratio_pco = pcolor(obs_L_grid_bound,obs_PA_grid_bound,hete_cps_ratio_plot);
set(hete_ratio_pco,'EdgeColor','none');
ylim(yrange); % xlim([2.4,5]);
colormap(colormap_viridis);
xlim(xrange);
clim(clim_range);
% clim_now = clim;
c = colorbar('FontSize',14);
c.Label.FontSize = 16;
c.LineWidth = 1.5;
c.Label.String = 'Cassini CPS Ratio';
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize', 14, 'Layer','top');
% xlabel('L-shell','FontSize',16);
ylabel('PA (deg)','FontSize',16);
set(gca,'TickDir','out');
set(gca,'YTick',[10:20:90]);

nexttile
pco_x = repmat(L_grid_edge_full,length(K_grid_edge),1);
pco_y = PA_dis';
pco_y(:,end) = pco_y(:,end-1);
pco_v = flux_dis_ratio';
pco = pcolor(pco_x,pco_y,pco_v); %mxn pa_grid_edge x L_grid_edge
set(pco,'EdgeColor','none');
colormap(colormap_viridis);
c = colorbar('FontSize',14);
c.Label.String = 'Simulation Flux Ratio';
c.Label.FontSize = 16;
c.LineWidth = 1.5;
clim(clim_range);
% clim([0.2,max(pco_v,[],'all')]);
ylim(yrange);
xlim(xrange);
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize', 14,'Layer','top');
% set(gca,'ColorScale','log');
% clim([1e-4,1e2]);
% xlabel('L-shell','FontSize',16);
ylabel('PA (deg)','FontSize',16);
set(gca,'TickDir','out');
set(gca,'YTick',[10:20:90]);

% cps or flux line
nexttile
errorbar(obs_L_grid,hete_cps_now_line,hete_cps_sem_line,'k-','LineWidth',2); hold on;
set(gca,'YScale','log');
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize', 14);
xlim([2.4,5.1]);
ylim([300,1.5e4]);
ylabel('CPS (#/s)','FontSize',16);
set(gca,'TickDir','out');

yyaxis right
enlarge = 200;
plot(L_grid_full,flux_dis_line.*enlarge,'b-','LineWidth',2);
set(gca,'YScale','log','YColor','b');
xlim(xrange);
ylim([300,1.5e4]);
ylabel(['Flux x ',num2str(enlarge)],'FontSize',16);
xlabel('L-shell','FontSize',16);

% window control
% screen size
set(0, 'Units', 'pixels');
screenSize = get(0, 'ScreenSize');
screenHeight = screenSize(4);

% set the height of figure to the screen height
fig.Units = 'pixels';
fig.Position(2) = 0;
fig.Position(4) = screenHeight; %575;
fig.Position(3) = 500; %670;
fig.Position(1) = (screenSize(3) - fig.Position(3)) / 2;

