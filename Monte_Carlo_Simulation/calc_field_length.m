% field line length
function res = calc_field_length(L,s)
res = L.*(asinh(sqrt(3).*s)./(2*sqrt(3))+s.*sqrt(1+3*s.^2)./2); %[Rs]
end