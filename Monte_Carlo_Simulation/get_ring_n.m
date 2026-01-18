% ring distribution
function res = get_ring_n(R)
res = zeros(1,size(R,2));
rho = sqrt(R(1,:).^2+R(2,:).^2);
z = R(3,:); %[Rs]
z_offset = -0.038;
rho_c = 3.95;
% radial
p_in = 16;
p_out = -10;
% vertical
w_in = [0.034,-0.066];
w_out = [0.034,0.044];
% offset
z0_param = [-0.005,0.018];

% inner
pickid = rho<=rho_c;
w = w_in(1)+w_in(2).*(rho(pickid)-rho_c);
Z0 = z0_param(1)+z0_param(2).*(rho(pickid)-rho_c)+z_offset;
res(pickid) = (rho(pickid)./rho_c).^p_in.*(w.^2./((z(pickid)-Z0).^2+w.^2));
% outer
pickid = rho>rho_c;
w = w_out(1)+w_out(2).*(rho(pickid)-rho_c);
Z0 = z0_param(1)+z0_param(2).*(rho(pickid)-rho_c)+z_offset;
res(pickid) = (rho(pickid)./rho_c).^p_out.*(w.^2./((z(pickid)-Z0).^2+w.^2));
end
