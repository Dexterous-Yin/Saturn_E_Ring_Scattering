% derive K from PA
function aeq = derive_K(PA,L,mapK,mapL,mapPA)
aeq = griddata(mapL,mapPA,mapK,L,PA);
% aeq = interp2(mapL,mapK,mapPA,L,Kq);
% Lid = find(L==mapL(:,1));
% PA = mapPA(Lid,:);
% K = mapK(Lid,:);
% aeq = interp1(K,PA,Kq); % [deg]
end