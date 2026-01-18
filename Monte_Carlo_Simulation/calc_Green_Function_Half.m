% set up the (L,K,mu) grids and calculate the corresponding physical quantities
% using Ye Model

% basic parameters
load('./Kmap.mat');
[params, ~, ~] = set_parameters('e');
ringparams = set_ring_parameters();
delL = 0.1;
L_grid = 2.4:delL:5;
L_grid_edge = min(L_grid)-delL/2:delL:max(L_grid)+delL/2;
% loss cone
Beq = params.Bs./L_grid.^3;
h_loss = params.Rs; lat_loss = acos(sqrt(h_loss./L_grid./params.Rs));
Bloss = params.Bs./L_grid.^3.*sqrt(1+3*sin(lat_loss).^2)./cos(lat_loss).^6;
lc_eq_grid = asin(sqrt(Beq./Bloss));

%% test K
PA_test = [3:1:90];
PA_test_edge = [3,(PA_test(1:end-1)+PA_test(2:end))/2,90];
K_test = zeros(1,length(PA_test));
for j = 1:length(PA_test)
    tempK = calc_K(PA_test(j),5,params);
    K_test(j) = tempK;
end
K_test_edge = zeros(1,length(PA_test_edge));
for j = 1:length(PA_test_edge)
    tempK = calc_K(PA_test_edge(j),5,params);
    K_test_edge(j) = tempK;
end
K_test = K_test(end:-1:1);
K_test_edge = K_test_edge(end:-1:1);

%% mu, K grids
mu_grid_min = 3.6;
mu_grid_max = 2.7e5;
mu_grid_num = 70; % 70 grid points are set, but only 61 are used, see del_grid next.
mu_grid = exp(interp1([1,mu_grid_num],[log(mu_grid_min),log(mu_grid_max)],1:1:mu_grid_num,'linear','extrap'));
mu_grid_edge = exp(interp1([1,mu_grid_num],[log(mu_grid_min),log(mu_grid_max)],0.5:1:mu_grid_num+0.5,'linear','extrap'));
K_grid = K_test; 
K_grid_edge = K_test_edge; 

%% corresponding E and PA grids
E_grid = zeros(length(L_grid),length(K_grid),length(mu_grid));
PA_grid = zeros(length(L_grid),length(K_grid),length(mu_grid));
for Li=1:length(L_grid)
    tempa = derive_PA(K_grid,L_grid(Li)+zeros(1,length(K_grid)),mapK,mapL,mapPA);
    tempaa = repmat(tempa',1,length(mu_grid));
    PA_grid(Li,:,:) = tempaa;
    disp(Li);
end
for Li=1:length(L_grid)
    for mui=1:length(mu_grid)
        tempE = derive_W(mu_grid(mui),squeeze(PA_grid(Li,:,mui)),L_grid(Li),params);
        E_grid(Li,:,mui) = tempE;
    end
    disp(Li);
end
E_grid_ori = E_grid;
PA_grid_ori = PA_grid;

E_grid_edge = zeros(length(L_grid),length(K_grid_edge),length(mu_grid_edge));
PA_grid_edge = zeros(length(L_grid),length(K_grid_edge),length(mu_grid_edge));
for Li=1:length(L_grid)
    tempa = derive_PA(K_grid_edge,L_grid(Li)+zeros(1,length(K_grid_edge)),mapK,mapL,mapPA);
    tempaa = repmat(tempa',1,length(mu_grid_edge));
    PA_grid_edge(Li,:,:) = tempaa;
    disp(Li);
end
for Li=1:length(L_grid)
    for mui=1:length(mu_grid_edge)
        tempE = derive_W(mu_grid_edge(mui),squeeze(PA_grid_edge(Li,:,mui)),L_grid(Li),params);
        E_grid_edge(Li,:,mui) = tempE;
    end
    disp(Li);
end

%% del grid
del_grid = zeros(length(L_grid),length(K_grid),length(mu_grid));
for Li=1:length(L_grid)
    delid = find(squeeze(PA_grid(Li,:,:))<lc_eq_grid(Li)/pi*180);
    del_grid(Li,delid) = 1;
end
for Ki=1:length(K_grid)
    for mui=1:length(mu_grid)
        temp = find(E_grid(:,Ki,mui)>=500 & E_grid(:,Ki,mui)<=10e3); % limit in 10MeV.
        if isempty(temp)
            del_grid(:,Ki,mui) = 1;
        end
    end
end

%% calculate collision number and trajectory length for particles at each L, PA
PA_set = 3:0.005:90;
coll_num = zeros(length(L_grid),length(PA_set),length(ringparams.size));
coll_length = zeros(length(L_grid),length(PA_set));
parfor Li=1:length(L_grid)
    tic
    [tempcoll_num,tempcoll_length] = calc_coll(L_grid(Li),PA_set,params,ringparams);
    coll_num(Li,:,:) = tempcoll_num';
    coll_length(Li,:) = tempcoll_length;
    toc
    disp(Li);
end

save('./Green_Function_First_Half.mat');
