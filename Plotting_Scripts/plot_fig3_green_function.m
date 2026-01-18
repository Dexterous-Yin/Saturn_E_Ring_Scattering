%{
This script is used to plot the Green function.
%}
% load matplotlib.pyplot
pyplt = py.importlib.import_module('matplotlib.pyplot');
cmap = pyplt.get_cmap('viridis');
cmap_vals = cmap(linspace(0, 1, 256));
cmap_mat = double(py.array.array('d', py.numpy.nditer(cmap_vals)));
cmap_mat = reshape(cmap_mat, [4, 256])';
colormapData = cmap_mat(:,1:3);
colormap_viridis = colormapData;

%% load basic settings
load('../Monte_Carlo_Simulation/Data/Green_Function_First_Half.mat');
% bounce period for each grid
[K_grid_now,L_grid_now,mu_grid_now] = meshgrid(K_grid,L_grid,mu_grid);
W_grid = E_grid.*params.keV;
gama_grid = W_grid./params.me/params.csq+1;
V_grid = sqrt((gama_grid.^2-1)*params.csq./gama_grid.^2);
Tb_grid = zeros(size(E_grid));
for Li = [12,17]
    for Ki = 1:length(K_grid)
        for mui = 1:length(mu_grid)
            coll_length_temp=interp1(PA_set,squeeze(coll_length(Li,:)),PA_grid(Li,Ki,mui));
            Tb_grid(Li,Ki,mui) = coll_length_temp/V_grid(Li,Ki,mui);
        end
    end
    disp(Li);
end

%% Green functions
dataprefix = '../Monte_Carlo_Simulation/Data/';
load([dataprefix,'Green_full.mat']);
n_particles = 1e6;
% modify for delted grids
for Li = [12,17]
    for Ki = 1:length(K_grid) % Ki in
        for mui = 1:length(mu_grid) % mui in
            if del_grid(Li,Ki,mui) == 1
                Green_full(Li,Ki,mui,Ki,mui) = n_particles;
            end
        end
    end
end
% calculate dedomega
dEdomega = zeros(size(E_grid)); % [keV*sr]
dEdsineq = zeros(size(E_grid)); % [keV*sin^2]
for Li=[12,17]
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

% modify the Green function for phase space density
Green_psd = zeros(size(Green_full));
counts_input = zeros(size(E_grid))+n_particles;
psd_input = counts_input./dEdsineq./p_sq_grid./Tb_grid;
for Ki = 1:length(K_grid) % Ki in
    for mui = 1:length(mu_grid) % mui in
        green_now = squeeze(Green_full(:,Ki,mui,:,:));
        green_now_psd = green_now./dEdsineq./p_sq_grid./Tb_grid;
        for Li = [12,17] %1:length(L_grid)
            green_ratio = squeeze(green_now_psd(Li,:,:))./psd_input(Li,Ki,mui);
            Green_psd(Li,Ki,mui,:,:) = green_ratio;
        end
    end
end
%%
fig = figure('Color',[1 1 1]);
layout = tiledlayout(3,2, 'TileSpacing', 'compact', 'Padding', 'compact');

% model
nexttile;
L = 2.4:0.01:5;
n = zeros(size(L));
n0 = 0.035;
for Li = 1:length(L)
    lat = -40:0.001:40;
    x = L(Li).*cosd(lat).^3;
    y = zeros(size(lat));
    z = L(Li).*cosd(lat).^2.*sind(lat);
    R = [x;y;z];
    n_temp = n0*get_ring_n(R);
    n(Li) = max(n_temp);
end
plot(L,n,'b-','LineWidth',2.5); hold on;
ylim([0,0.04]); xlim([2.4,5]);
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize', 14);
set(gca,'TickDir','out','TickLength',[0.015,0.015],'YTick',[0:0.01:0.04]);
xlabel('L-shell','FontSize',16);
ylabel('n (m^{-3})','FontSize',16);
plot([4,4],[0,0.05],'k--','LineWidth',1.5);
plot([3.5,3.5],[0,0.05],'k--','LineWidth',1.5);

%
% collusion frequency
nexttile
time_norm = 1;
W_temp = 1e3.*params.keV;
gama_temp = W_temp./params.me/params.csq+1;
v_temp = sqrt((gama_temp.^2-1)*params.csq./gama_temp.^2);
coll_time_temp = coll_length./v_temp;
for tempL = [12,17]
    coll_num_plot = sum(squeeze(coll_num(tempL,:,:)),2).*(time_norm./coll_time_temp(tempL,:))';
    if tempL == 12
        plot(PA_set,coll_num_plot,'b-','LineWidth',2.5); hold on;
    else
        plot(PA_set,coll_num_plot,'-','Color','#D95319','LineWidth',2.5); hold on;
    end
end
xlim([3,90]);
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize', 14);
set(gca,'TickDir','out','TickLength',[0.015,0.015]);
xlabel('PA (deg)','FontSize',16);
ylabel('Collion Frequency (s^{-1})','FontSize',16);
set(gca,'YScale','log');
ylim([2e-7,5e-5]);

% window control
% screen size
set(0, 'Units', 'pixels');
screenSize = get(0, 'ScreenSize');
screenHeight = screenSize(4);
screenWidth = screenSize(3);

% set the height of figure to the screen height
fig.Units = 'pixels';
fig.Position(4) = screenHeight;
fig.Position(3) = 680; %screenWidth;
fig.Position(1) = (screenSize(3) - fig.Position(3)) / 2;
fig.Position(2) = (screenHeight - fig.Position(4)) / 2;

% nexttile
plot_params = [12,6,33; ...
    17,6,33; ...
    12,53,27; ...
    17,53,27];

for ploti = 1:length(plot_params(:,1))
    Li = plot_params(ploti,1);
    Ki = plot_params(ploti,2);
    mui = plot_params(ploti,3);

    xlim_range = [70,2.5e3];

    tempGreen = squeeze(Green_psd(Li,Ki,mui,:,:));
    pick_PA_grid_edge = squeeze(PA_grid_edge(Li,:,:));
    pick_E_grid_edge = squeeze(E_grid_edge(Li,:,:));
    pick_Green_grid = tempGreen;

    pick_E_grid = squeeze(E_grid(Li,:,:));
    pick_PA_grid = squeeze(PA_grid(Li,:,:));

    nexttile;
    % plot_PA_grid = pick_PA_grid(end:-1:1,1)';
    plot_Green_grid = pick_Green_grid(end:-1:1,:);
    plot_Green_grid_full = zeros(length(pick_PA_grid_edge),length(mu_grid_edge));
    plot_Green_grid_full(1:length(pick_PA_grid_edge)-1,1:length(mu_grid_edge)-1) = plot_Green_grid;
    % pco1 = pcolor(mu_grid,plot_PA_grid,plot_Green_grid);
    pco1 = pcolor(mu_grid_edge,pick_PA_grid_edge(end:-1:1,1)',plot_Green_grid_full);
    hold on;
    % plot(mu_grid(mui),K_grid(Ki),'LineStyle','-','LineWidth',4,'Marker','x','MarkerSize',14,'MarkerFaceColor','m','MarkerEdgeColor','m');
    set(pco1,'EdgeColor','none');
    set(gca,'XScale','log','ColorScale','log');
    % set(gca,'YScale','log');
    % ylim([1e-4,0.2]);
    if ploti<=2
        ylim([70,90]);
    else
        ylim([30,50]);
    end
    ylim_now = ylim;
    colormap(gca,colormap_viridis);
    % clim([5e-7,max(pick_Green_grid,[],'all','omitnan')]);
    clim([1e-6,1]);
    clim_now = clim;
    c = colorbar('FontSize',14);
    c.Label.FontSize = 16;
    c.LineWidth = 1.5;
    c.TickLength = 0.015;
    % c.Ticks = [1e-5,1e-3,1e-1];
    % c.Label.String = 'PSD Ratio';
    set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize', 14,'Color',[0.8,0.8,0.8], 'Layer','top');
    set(gca,'TickDir','out','TickLength',[0.015,0.015]);
    xlabel('mu (MeV/G)','FontSize',16);
    ylabel('PA (deg)','FontSize',16);
    xlim(xlim_range);
    xticks = 1:1:4;
    xticklabels = {'10^1','10^2','10^3',''};
    set(gca,'XTick',10.^xticks,'XTickLabel',xticklabels);

    % add equal E lines
    line_E = [E_grid(Li,Ki,mui)];
    for li = 1:length(line_E)
        line_Enow = line_E(li);
        line_pa = ylim_now(1):0.1:ylim_now(2);
        line_mu = calc_mu(line_Enow,line_pa,L_grid(Li),params);
        plot(line_mu,line_pa,'w--','LineWidth',1);
        text(line_mu(end),line_pa(end)-2,[num2str(E_grid(Li,Ki,mui)/1e3,'%.2f'),' MeV'],'Color','w','FontSize',16);
        text(line_mu(1),line_pa(1)+2,[num2str(tempGreen(Ki,mui),'%.3f')],'Color','w','FontSize',16);
    end

    hold on;
    plot(mu_grid(mui),pick_PA_grid(Ki,1),'LineStyle','none','Marker','o','MarkerSize',8,'MarkerFaceColor',"#cccccc",'MarkerEdgeColor',"#cccccc");


end