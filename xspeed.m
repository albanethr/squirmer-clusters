function [ux] = xspeed(h,D,beta,theta)
% horizontal translation speed (x-direction) for two squirmers of strength beta at height h above a no slip wall and distance D from one another 
% The translation speed used is the speed of the surrounding fluid evaluated at the centre of the sphere

R = 1;
s = 1/2*R^3; 
p = -3/4*beta*R;

u1 = cos(theta); 

y = 2.*h; 
x = -D; 
r = sqrt(x.^2+y.^2);

ufd2i = p .*(x./(r.^3) + ( (-3.*x.^3 + 18.*h.^2 .*x - 18.*h.*y.*x )./(r.^5) + (30.*h.*x.^3.*(y-h) )./(r.^7) ).*cos(theta).^2+ ((-6.*h.*x.*y -3.*x.*y.^2+6.*h.^2.*x)./(r.^5) + (30.*h.*x.*y.^2.*(y- h))./(r.^7) ).* sin(theta).^2 - (( 6.*x.^2.*y-12.*h.^2.*y + 12.*h.*y.^2)./(r.^5) + ( 60.*h.*x.^2.*y.*(h-y) )./(r.^7) ).* sin(theta) .* cos(theta) );
 usd2i = s.*(-cos(theta).*( r.^-3 + (-3.*x.^2-6.*y.^2+6.*h.*y)./(r.^5) + (30.*h.*x.^2.*y.*(y-h))./(r.^7))+ sin(theta).*((9.*x.*y-6.*h.*x)./(r.^5) + (30.*x.*y.^2.*(h-y))./(r.^7)));


y = 0; 
r = sqrt(x.^2+y.^2);

usd2=  s./r.^3 .*( cos(theta) + 3.*x./r.^2.*(-x.*cos(theta)+y.*sin(theta)));
ufd2 = - p.*x./(r.^3) .*(1 - 3./(r.^2) .* (x.^2.* cos(theta).^2 + y.^2.*sin(theta).^2) ) ;

usd1i = - s.*cos(theta)./(4.*h.^3); 
ufd1i =  3.*p .*sin(2.*theta)./(8.*h.^2);

ux =  (u1 + ufd1i + usd1i + ufd2i +  ufd2  + usd2 + usd2i); 

urep = - ux*(0.5./abs(D-2*R) );

ux = ux + urep; 

end