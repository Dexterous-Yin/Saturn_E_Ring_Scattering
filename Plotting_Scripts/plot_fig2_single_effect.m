%{
This script is used to plot the Geant4 results for the effect of a single ice grain.
%}

% load matplotlib.pyplot
pyplt = py.importlib.import_module('matplotlib.pyplot');
cmap = pyplt.get_cmap('viridis');
cmap_vals = cmap(linspace(0, 1, 256));
cmap_mat = double(py.array.array('d', py.numpy.nditer(cmap_vals)));
cmap_mat = reshape(cmap_mat, [4, 256])';
colormapData = cmap_mat(:,1:3);
colormap_turbo = colormapData; %turbo;

%% PART1: plot probability distributions

% load Geant4 results
load('../Geant4_Simulation/full_Geant4_simulation.mat');

% plot
fig = figure('Color',[1 1 1]);
layout = tiledlayout(3,2, 'TileSpacing', 'tight', 'Padding', 'compact');

% effect of a single particle on electrons with different energies
% delta E
nexttile(1,[1,1]);
pco_delE = pcolor(delE_edge,G4E_set_edge,delE_plot);
set(pco_delE,'EdgeColor','none');
set(gca,'YScale','log');
ylim([200,5e3]); xlim([-3,0]);
clim([0,0.03]);
set(gca,'Colormap',colormap_turbo);
cbar = colorbar('FontSize', 16);
cbar.Label.FontSize = 16;
cbar.LineWidth = 1.5;
cbar.Label.String = 'Probility';
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize', 16, 'Layer',' top');
set(gca,'YTick',[250,5e2,1e3,2e3,4e3]);
xlabel('Delta E (keV)','FontSize',16); 
ylabel({'Ein (keV)'},'FontSize',16);
title([num2str(ringparams.size(20),'%.0f'),'-um ice'],'FontSize',16);
set(gca,'TickDir','out');

% relative delta E
nexttile(3,[1,1]);
pco_relE = pcolor(relE_edge*100,G4E_set_edge,relE_plot);
set(pco_relE,'EdgeColor','none');
set(gca,'YScale','log');
ylim([200,5e3]); xlim([-3e-3,0]*100);
clim([1e-3,0.15]);
set(gca,'ColorScale','log');
set(gca,'Colormap',colormap_turbo);
cbar = colorbar('FontSize', 16);
cbar.Label.FontSize = 16;
cbar.LineWidth = 1.5;
cbar.Label.String = 'Probility';
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize', 16, 'Layer',' top');
set(gca,'YTick',[250,5e2,1e3,2e3,4e3]);
xlabel('Relative dE (%)','FontSize',16); 
ylabel({'Ein (keV)'},'FontSize',16);
set(gca,'TickDir','out');

% relative delta E
nexttile(5,[1,1]);
pco_delPA = pcolor(delPA_edge,G4E_set_edge,delPA_plot);
set(pco_delPA,'EdgeColor','none');
set(gca,'YScale','log');
ylim([200,5e3]); xlim([0,10]);
clim([1e-3,0.2]);
set(gca,'ColorScale','log');
set(gca,'Colormap',colormap_turbo);
cbar = colorbar('FontSize', 16);
cbar.Label.FontSize = 16;
cbar.LineWidth = 1.5;
cbar.Label.String = 'Probility';
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize', 16, 'Layer',' top');
set(gca,'YTick',[250,5e2,1e3,2e3,4e3]);
xlabel('Delta Angle (deg)','FontSize',16); 
ylabel({'Ein (keV)'},'FontSize',18);
set(gca,'TickDir','out');
set(gca,'XTick',[0:2:10]);

% window control
% screen size
set(0, 'Units', 'pixels');
screenSize = get(0, 'ScreenSize');
screenWidth = screenSize(3);
screenHeight = screenSize(4);

% set the height of figure to the screen height
fig.Units = 'pixels'; 
fig.Position(2) = 0; 
fig.Position(4) = screenHeight;
fig.Position(3) = 800; %screenWidth;
fig.Position(1) = (screenSize(3) - fig.Position(3)) / 2;

%
% effect of particle with different sizes on electrons
% delta E
nexttile(2,[1,1]);
pco_delE = pcolor(fixE_delE_edge,ringparams.size_bound,fixE_delE_plot);
set(pco_delE,'EdgeColor','none');
set(gca,'YScale','log');
ylim([1,10]); xlim([-3,0]);
% clim([0,0.06]);
clim([1e-3,0.06]); 
set(gca,'ColorScale','log');
set(gca,'Colormap',colormap_turbo);
cbar = colorbar('FontSize', 16);
cbar.Label.FontSize = 16;
cbar.LineWidth = 1.5;
cbar.Label.String = 'Probility';
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize', 16, 'Layer',' top');
set(gca,'YTick',[1,2,4,8]); %set(gca,'YTick',[1,2,4,8]);
xlabel('Delta E (keV)','FontSize',16); 
ylabel({'Size (um)'},'FontSize',16);
title(['1-MeV e-'],'FontSize',16);
set(gca,'TickDir','out');
%
% relative delta E
nexttile(4,[1,1]);
pco_relE = pcolor(fixE_relE_edge*100,ringparams.size_bound,fixE_relE_plot);
set(pco_relE,'EdgeColor','none');
set(gca,'YScale','log');
ylim([1,10]); xlim([-3e-3,0]*100);
% clim([0,0.04]);
clim([1e-3,0.06]); 
set(gca,'ColorScale','log');
set(gca,'Colormap',colormap_turbo);
cbar = colorbar('FontSize', 16);
cbar.Label.FontSize = 16;
cbar.LineWidth = 1.5;
cbar.Label.String = 'Probility';
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize', 16, 'Layer',' top');
set(gca,'YTick',[1,2,4,8]);
xlabel('Relative dE (%)','FontSize',16); 
ylabel({'Size (um)'},'FontSize',16);
set(gca,'TickDir','out');

% relative delta E
nexttile(6,[1,1]);
pco_delPA = pcolor(fixE_delPA_edge,ringparams.size_bound,fixE_delPA_plot);
set(pco_delPA,'EdgeColor','none');
set(gca,'YScale','log');
ylim([1,10]); xlim([0,10]);
% clim([0,0.1]);
clim([1e-3,0.1]); 
set(gca,'ColorScale','log');
set(gca,'Colormap',colormap_turbo);
cbar = colorbar('FontSize', 16);
cbar.Label.FontSize = 16;
cbar.LineWidth = 1.5;
cbar.Label.String = 'Probility';
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize', 16, 'Layer',' top');
set(gca,'YTick',[1,2,4,8]);
xlabel('Delta Angle (deg)','FontSize',16); 
ylabel({'Size (um)'},'FontSize',16);
set(gca,'TickDir','out');
set(gca,'XTick',[0:2:10]);

%% PART2 (bottom panels): Results with physical processes
% 1 MeV 3um
filename = '../Geant4_Simulation/sample_output/Record_1000.0000_6.00.mat'; 
load(filename);
E_process_1 = E;
PA_process_1 = PA;
HasEBrem_process_1 = HasEBrem;
HasEIoni_process_1 = HasEIoni;
clearvars E PA HasEBrem HasEIoni
% 1 MeV 10um
filename = '../Geant4_Simulation/sample_output/Record_1000.0000_20.00.mat'; 
load(filename);
E_process_2 = E;
PA_process_2 = PA;
HasEBrem_process_2 = HasEBrem;
HasEIoni_process_2 = HasEIoni;
clearvars E PA HasEBrem HasEIoni

% plot
fig = figure('Color',[1 1 1]);
layout = tiledlayout(1,2, 'TileSpacing', 'tight', 'Padding', 'compact');

nexttile
bothid = find(HasEBrem_process_1 == 1 & HasEIoni_process_1 == 1);
ebremid = find(HasEBrem_process_1 == 1 & HasEIoni_process_1 ==0);
eioniid = find(HasEBrem_process_1 == 0 & HasEIoni_process_1 == 1);
neitherid = find(HasEBrem_process_1 == 0 & HasEIoni_process_1 == 0);
% plot
dot_size = 40;
scatter(E_process_1(neitherid),PA_process_1(neitherid),dot_size,'filled','MarkerFaceColor','#7f7f7f','MarkerFaceAlpha',0.7,'DisplayName','Neither'); hold on;
scatter(E_process_1(ebremid),PA_process_1(ebremid),dot_size,'b','filled','MarkerFaceAlpha',0.7,'DisplayName','eBrem'); hold on;
scatter(E_process_1(eioniid),PA_process_1(eioniid),dot_size,'filled','MarkerFaceColor','#008b4c','MarkerFaceAlpha',0.7,'DisplayName','eIoni');
box on;
% legend show;
set(gca,'XMinorTick','on','YMinorTick','off','LineWidth',2,'FontSize', 16);
xlim_now = xlim;
xlim([450,1e3+10]);
ylim([-3,185]); set(gca,'YTick',[0:20:180]);
xlabel('dE (keV)','FontSize',16);
ylabel('Del Angle (deg)','FontSize',16);
set(gca,'TickDir','out');
title('1-MeV e- & 3-um ice','FontSize',16);
% set(gca,'YTickLabel',{});
% set(gca,'XTickLabel',{});
% set(gca,'XTick',[],'YTick',[],'Box','off','LineWidth',0.5); 


nexttile
bothid = find(HasEBrem_process_2 == 1 & HasEIoni_process_2 == 1);
ebremid = find(HasEBrem_process_2 == 1 & HasEIoni_process_2 ==0);
eioniid = find(HasEBrem_process_2 == 0 & HasEIoni_process_2 == 1);
neitherid = find(HasEBrem_process_2 == 0 & HasEIoni_process_2 == 0);
% plot
dot_size = 40;
scatter(E_process_2(neitherid),PA_process_2(neitherid),dot_size,'filled','MarkerFaceColor','#7f7f7f','MarkerFaceAlpha',0.7,'DisplayName','Neither'); hold on;
scatter(E_process_2(ebremid),PA_process_2(ebremid),dot_size,'b','filled','MarkerFaceAlpha',0.7,'DisplayName','eBrem'); hold on;
scatter(E_process_2(eioniid),PA_process_2(eioniid),dot_size,'filled','MarkerFaceColor','#008b4c','MarkerFaceAlpha',0.7,'DisplayName','eIoni');
box on;
% legend show;
set(gca,'XMinorTick','on','YMinorTick','off','LineWidth',2,'FontSize', 16);
xlim_now = xlim;
xlim([450,1e3+10]);
ylim([-3,185]); set(gca,'YTick',[0:20:180]);
xlabel('dE (keV)','FontSize',16);
ylabel('Del Angle (deg)','FontSize',16);
set(gca,'TickDir','out');
title('1-MeV e- & 10-um ice','FontSize',16);
% set(gca,'YTickLabel',{});
% set(gca,'XTickLabel',{});
% set(gca,'XTick',[],'YTick',[],'Box','off','LineWidth',0.5); 

% window control    
% screen size
set(0, 'Units', 'pixels');
screenSize = get(0, 'ScreenSize');
screenWidth = screenSize(3);
screenHeight = screenSize(4);

% set the height of figure to the screen height
fig.Units = 'pixels'; 
fig.Position(2) = screenHeight*1/4; 
fig.Position(4) = screenHeight*1/4;
fig.Position(3) = 800;
fig.Position(1) = (screenSize(3) - fig.Position(3)) / 2;