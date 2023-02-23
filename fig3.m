function [] = fig3()
% Returns the equlibrium position (height and orientation) of a single
% squirmer above a boundary if they exist
% The figure is the phase diagram in fig 3 

figure
hold on;

% range of explored \alpha and \beta
nalpha = 1000; nbeta = 2000; 
alpha = linspace(0.3,0.999,nalpha); B = linspace(-2,2,nbeta);
[X,Y] = meshgrid(B,alpha); height = 0.*X -1;

for i=1:nbeta
    for j = 1:nalpha
        beta = B(i);
        v0 = alpha(j); vg =1; 

        % in the case of pullers 
        if beta > 0 
            p = [-(v0+vg)/2 ; -9*v0*beta/16; 9/8*vg; v0-vg];
            hh = roots(p).^-1;
            if isreal(max(hh))
                if max(hh)>0 && polyval(polyder(p),max(hh).^-1) > 0, height(j,i)= max(hh); end
            end
         
        % pushers    
        else 
            % upright pusher case 
            vg =1; 
            p = [-(v0+vg)/2 ; -9*v0*beta/16; 9/8*vg; v0-vg];
            hh = roots(p).^-1;
            if isreal(max(hh))
                if max(hh)>0 && polyval(polyder(p),max(hh).^-1) > 0 && abs(beta) < 2/3 /max(hh) 
                    height(j,i)= max(hh); 
                end
            end
            % find the roots of the polynom in eq (1.9) 
            p = 1e6.*[ -1/(24*beta); - 1/(2*v0); 9*beta/32 ; 9/(8 * v0) - 2/ (3*beta); -1/v0 ];
            hse = roots(p).^-1;
            for k=1:4
                if isreal(hse(k)) % check if it has real roots
                    if hse(k)> 0 && polyval(polyder(p),hse(k).^-1) > 0 && polyval(p,hse(k).^-1) < 1e-6
                        if abs(beta) > 2/3 /hse(k) % condition for tilted stability 
                        height(j,i) = hse(k);    
                        end
                    end
                end
            end
            
        end
    end
end



surf(X,Y,height,'EdgeColor','none');
zlim([1,max(max(height))+1]); % plot only the equilibria that exist
ylim([0.3,1]);

view(2) 

ylabel('\alpha')
xlabel('\beta')

caxis([1,7])

% analytical boundary of the tilted equilibrium region 
beta = linspace(-2/3,0);
plot3(beta, -4*(27*beta.^3-27.*beta-16)./(27*beta.^3+64), 100+0.*beta,'color',[0.64,0.08,0.18],'LineWidth',1.2);
beta = linspace(-1.5,-2/3);
plot3(beta,36.*beta./(-68+27.*beta.^2),2+ 0.*beta,'color',[0.64,0.08,0.18],'LineWidth',1.2); 

% h = 2 line for tilted
beta = linspace(-1,-0.345);
plot3(beta,22.*beta./(-16+9.*beta.^2),20+ 0.*beta,'--','color',[0.64,0.08,0.18],'LineWidth',1); 

% h = 2 line for upright 
beta = linspace(-1/3,2); 
plot3(beta,-(32./(3.* (-20 + 3 .*beta))),20+ 0.*beta,'--','color',[0.64,0.08,0.18],'LineWidth',1); 

% boundaries for the lubrication analysis 
plot3(-1 + 0.* linspace(0.3,1), linspace(0.3,0.9), 1+0.* linspace(0.3,1),'--k','LineWidth',1); 
plot3(1 + 0.* linspace(0.3,1), linspace(0.3,0.62), 1.00001+ 0.* linspace(0.3,1),'--k','LineWidth',1); 

% pusher/puller distinction and experimental alpha
plot3(0 + 0.* linspace(0.3,1), linspace(0.3,1), max(max(height))+1 + 0.* linspace(0.3,1),':k','LineWidth',1); 
plot3(linspace(-2,2), .77+0*linspace(0.3,1), max(max(height))+1 + 0.* linspace(0.3,1),':k','LineWidth',1); 
end



