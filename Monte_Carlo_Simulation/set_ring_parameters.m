function ringparams = set_ring_parameters()
    % ring parameters
    size_start = 1; %[um]
    size_end = 10;
    size_num = 40;
    srange = exp(interp1([1,size_num],[log(size_start),log(size_end)],1:1:size_num));
    srange_bound = exp(interp1([1,size_num],[log(size_start),log(size_end)],0.5:1:size_num+0.5,'linear','extrap'));
    ringparams.size = srange;
    ringparams.size_bound = srange_bound;
    q = -4; % coefficient for differential size distribution
    ntot = 0.0321; %0.1832;
    ringparams.ringn = ntot*(srange_bound(1:end-1).^(q+1)-srange_bound(2:end).^(q+1));
    ringparams.cross = pi*(ringparams.size*1e-6).^2;
end

