% new PA
function y=calc_newPA(dtheta,theta,phi)
% calculate the new pitch angle
% Note the scattered angle (dtheta) is relative to the direction of
% initial momentum, not the magnetic field.
% theta: initial pitch angle
% dtheta: scattered angle
% phi: random azimuthal scatter angle
y=acos(cos(theta).*cos(dtheta)+sin(theta).*sin(dtheta).*cos(phi));
end