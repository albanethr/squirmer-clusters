function [uy] = yspeed(h,D,beta,theta)
% vertical translation speed (y-direction) for two squirmers of strength beta at height h above a no slip wall and distance D from one another 
% The translation speed used is the speed of the surrounding fluid evaluated at the centre of the sphere

R = 1;

% strenght of the hydrodynamic singularities 
s = 1/2*R^3; 
p = -3/4*beta*R^2;

alpha = 20/26; %v0/vg

y = 2.*h; 
x = -D; 
r = sqrt(x.^2+y.^2);

u1 = sin(theta);

% decompose between the different components
usd1i = - s.*sin(theta)./(h.^3) ;
ufd1i = - 3.*p.*(1-3.*sin(theta).^2)./(8.*h.^2) ;

ufd2i = p* ((y-2.*h)./(r.^3) +((6.*h.*x.^2-3.*x.^2.*y +6.*h.^2.*y - 6.*h.*y.^2)./(r.^5) + (30.*h.*x.^2.*y.*(y-h ))./(r.^7)) .*cos(theta).^2   + ((-3.*y.^3 + 18.*h.^2.*y - 12.*h.*y.^2)./(r.^5) + (30.*h.*y.^3.*(y-h))./(r.^7) ).*sin(theta).^2 -  ((9.*x.*y.^2 -12.*h.^2.*x -3.*x.*y.^2)./(r.^5) + (60.*h.*y.^2.*x.*(h-y))./(r.^7) ).*sin(theta).*cos(theta));
usd2i = s *( -cos(theta).*(( -3.*x.*y + 6.*h.*x)./(r.^5) + ( 30.*x.*y.^2.*(y-h) )./(r.^7) ) + sin(theta).*( r.^-3 + (15.*y.^2 - 18.*h.*y )./(r.^5) + (30.*y.^3.*(h-y))./(r.^7) ));

y = 0; 
r = sqrt(x.^2+y.^2);

usd2=  s./r.^3 .*( -sin(theta) + 3.*y./r.^2.*(-x.*cos(theta)+y.*sin(theta)));
ufd2 = - p.*y./(r.^3) .*(1 - 3./(r.^2) .* (x.^2.* cos(theta).^2 + y.^2.*sin(theta).^2) ) ;

ug = - alpha^-1 * ( 1- 9/8.*(R./h) +1/2*(R./h).^3);

uy = (ug + ( u1 + ufd1i + usd1i)) +  ( usd2 + ufd2 + ufd2i + usd2i ) ;


end