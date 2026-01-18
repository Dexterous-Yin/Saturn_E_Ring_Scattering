% derive pitch angle from K
function aeq = derive_PA(Kq,L,mapK,mapL,mapPA)
aeq = griddata(mapL,mapK,mapPA,L,Kq);
% aeq = interp2(mapL,mapK,mapPA,L,Kq);
% Lid = find(L==mapL(:,1));
% PA = mapPA(Lid,:);
% K = mapK(Lid,:);
% aeq = interp1(K,PA,Kq); % [deg]
end