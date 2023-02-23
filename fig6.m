function [] = fig6()
% dynamics of two approaching squirmers starting from D at onesquirmer equilibrium
% uses the auxiliary function xspeed, yspeed, and rotationrate the
% translation and reorientation speed respectively
% returns the plot (b,c,f,g) of figure 6 for the chosen range of beta and
% initial separation D


% Chose range of beta 
bbeta =  linspace(0,2,8); % range of pullers 
%bbeta = linspace(-0.7,-0.2,7); % range of pushers 

% Chose the initial separation D 
D0 = 5; 

% other parameters 
R = 1;
vg = 26;
v0=20; vga = vg/20;


cm = colormap(hot(length(bbeta)+5));

for j = 1:length(bbeta)
    
    D = D0;
        
    figure(1);
    hold on;
    
    beta = bbeta(j);
    
    % find the starting single squirmer equilibrium
    h0=1; theta0=pi/2;
    
    if beta < -0.13
        h0 = fsolve(@(h) 1e6*(27* beta^2* h^2 - 4 *R* (-8 + 16* h^3 + 9* R) - 12* beta* h *(4 - 9* h^2 + 8* h^3)* vga)/(96* beta *h^4), 5);
        theta0 = - acos(-2*R/(3*beta*h0)) +pi/2;
    end
    
    if ~isreal(h0)|| ~isreal(theta0) || beta >= -0.13
        p = [-(v0+vg)/2 ; -9*v0*beta/16; 9/8*vg; v0-vg];
        htp = roots(p).^-1;
        for kj=1:length(htp)
            if isreal(htp(kj))
                if htp(kj) > 1.01 && polyval(polyder(p),htp(kj).^-1) > 0, h0=htp(kj);
                end
            end
        end
        theta0 = pi/2;
        disp('vertical');
    end
    disp(['h = ',num2str(h0),', theta = ',num2str(theta0)]);
    theta=theta0; h = h0; 
    
    dt = 0.01;
    Nt = 100000;
    
    noise = 0;
    
    tth = zeros(1,Nt); DD = zeros(1,Nt); hh=zeros(1,Nt);
    
    % Dynamics time-stepping 

    for i = 1:Nt
        
        thetan = theta + rotationrate(h,D,beta,theta)*dt + dt*noise*(rand()-1/2);
        Dn = D - 2 * dt * xspeed(h,D,beta,theta) + dt*noise*(rand()-1/2);
        hn = h + dt * yspeed(h,D,beta,theta) + dt*noise*(rand()-1/2);
        
        % collision 
        if D < 2
            break
        end
        
        theta = thetan; D = Dn; h = hn;
        
        if ~isempty(theta)
            tth(i) = theta; DD(i) = D; hh(i) = hn;
            % plot(D,h,'.r','Markersize',8,'Color', cm(i,:));
        end
        
    end
  
    subplot(2,1,1); hold on;
    hh = hh(DD>0);  DD = DD(DD>0);
    plot(DD,hh,'.-','Color', cm(j,:),'displayname',num2str(beta),'LineWidth',1);
    subplot(2,1,2); hold on; 
    plot(DD,tth(DD>0),'.-','Color', cm(j,:),'displayname',num2str(beta),'LineWidth',1);
      
end

end

