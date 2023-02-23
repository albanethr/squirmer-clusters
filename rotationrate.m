function [Om] = rotationrate(h,D,beta,theta)
% rotationrate for two squirmers at height h, distance D
% Om = 1/2 omega (the rotation of a spherical object is half that of the surrounding fluid) 

R = 1;
s = 1/2*R^3; 
p = -3/4*beta*R;

R = (D.^2+4.*h.^2).^.5;

Om = 6.*p.*h.*D./R.^7.*(-20.*h.*D.*sin(theta).*cos(theta)+ 2.*(D.^2-6.*h.^2).*cos(theta).^2 + (8.*h.^2-3.*D.^2).* sin(theta).^2)+ 6.*p.*h.*D.*R.^-5+ 3.*p.*sin(theta).*cos(theta).* (R.^-3+D.^-3 - h.^-3./8) + 6.*s.*R.^-7.* (D.*(D.^2-16.*h.^2).*sin(theta)+ 8.*h.*(D.^2 - h.^2).*cos(theta)) + 3.*s./(8.*h.^4).*cos(theta);

end

