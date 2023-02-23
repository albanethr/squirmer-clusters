function [] = partialcircle()
% This function explores the dynamics of a partial circle of squirmers in the lubrication regime 
% it is used in sec 3(g) 
% The figure displays the orientation of each swimmer as an arrow which goes from black to yellow in time 
% The simulation stops if the circle is unstable and displays the time at which it broke apart 

% set the properties of the circle 
N = 6;                  % number of squirmers
beta = -.2;             %squirmer parameter
h= 1.01;                % height above the wall (if h>2, this is the bulk configuration) 
wall2= 0;               % one for two walls 
epsis = 0.5;            % epsilon: separation between closest neighbours 
theta0 = 0;             % initial polar angle 


% Numerical parameters 
T = 5000;               % length of simulation
dt = 0.001;             % time step 
R = 1;
v0 = 1;
contnoise = 0.0001;     % noise at each time step
ininoise = 0;           % noise in the initial configuration 



X = zeros(3,N,T);
E = zeros(3,N,T);

cm = colormap(hot(floor(T*1.2)));

if N > 7, disp('N too large'); return; end

% circle configuration
D = 2+epsis; 
X(:,1:N-1,1) = [D*cos((1:1:N-1)*2*pi./6);h.*ones(1,(N-1));D*sin((1:1:(N-1))*2*pi./6)];
E(:,1:N-1,1) = [-cos(theta0).*cos((1:1:(N-1))*2*pi./6);sin(theta0).*ones(1,(N-1));-cos(theta0).*sin((1:1:(N-1))*2*pi./6)];
X(:,N,1)= [0;h;0];
E(:,N,1)= [0; 1 ; 0]; % rand(1,3);%[0;-1;0];
xlim([-D-1,D+1]); ylim([-D-1,D+1]); plot(X(1,1:N,1), X(3,1:N,1), '.','MarkerSize',12); daspect([1 1 1]); %return

% central orientation
et = normc(repmat([sum(X(1,:,1))./N;sum(X(3,:,1))./N],[1,N])-[X(1,:,1);X(3,:,1)]);
%plot(et(1),et(2),'.','MarkerSize',12)
E(:,:,1) = [cos(theta0).*et(1,:); sin(theta0).*ones(1,N); cos(theta0).*et(2,:)];

% noise
E(:,:,1) = E(:,:,1)+ [2*ininoise.*(rand(1,N)-1/2);0.*(rand(1,N)-1/2);2*ininoise.*(rand(1,N)-1/2)];
%E(:,:,1) =  [0.1,0.1;0,0;0.5,0.5];

E(:,:,1)= E(:,:,1)./repmat(vecnorm(E(:,:,1)),3,1);


figure(1)
clf
hold on;

quiver3(X(1,:,1),X(3,:,1),X(2,:,1),E(1,:,1),E(3,:,1),E(2,:,1),'LineWidth',1,'Color',cm(1,:));
drawnow();
xlabel('x'); ylabel('z');

unstable = 0; 

for t = 2:T
    for i = 1:N
        Omi = zeros(3,1);
        for j = [1:i-1, i+1:N]
            epsij = norm(X(:,i,t-1)-X(:,j,t-1)) - 2;
            if epsij < 1
                %disp([i,j])
                eparih = (X(:,j,t-1)-X(:,i,t-1))./norm(X(:,j,t-1)-X(:,i,t-1));

                %Omi1 = Omi + log(1./epsij)*12/(40*R)*v0*( (1+beta*dot(E(:,i,t-1), eparih) ) *cross( eparih, E(:,i,t-1)-dot(E(:,i,t-1), eparih)* eparih));
                Omi = Omi + log(1./epsij)*12/(40*R)*v0*( (1+beta*dot(E(:,i,t-1), eparih) ) *cross( eparih, E(:,i,t-1)));
                %Omi1 = Omi + log(1./epsij)*3/(40*R)*v0* ( (1+beta*dot(E(:,j,t-1), -eparih) ) *cross(-eparih, E(:,j,t-1)-dot(E(:,j,t-1),-eparih)*(-eparih)));
                Omi = Omi + log(1./epsij)*3/(40*R)*v0* ( (1+beta*dot(E(:,j,t-1), -eparih) ) *cross(-eparih, E(:,j,t-1)));

                if isnan(Omi), return; end
            end
        end
        E(:,i,t) = E(:,i,t-1)+dt.*cross(Omi,E(:,i,t-1));
        X(:,i,t) = X(:,i,t-1);
    end

        
    if h<2
        for i = 1:N
            Omw = log(1./(h-1))*24/(40*R)*v0* (1 - beta * E(2,i,t-1)) *cross([0;-1;0],[E(1,i,t-1);0;E(3,i,t-1)]);
            E(:,i,t) = E(:,i,t) + dt.*cross(Omw,E(:,i,t-1));
            if wall2==1
                Omw2 = log(1./(h-1))*24/(40*R)*v0* (1 + beta * E(2,i,t-1)) *cross([0;1;0],[E(1,i,t-1);0;E(3,i,t-1)]);
                E(:,i,t) = E(:,i,t) + dt.*cross(Omw2,E(:,i,t-1));
            end
        end
    end

    % 2D
    %E(2,:,t)= 0.*E(2,:,t);

    E(:,:,t) = E(:,:,t) + contnoise.*(2*rand(3,N)-1/2);

    E(:,:,t)= E(:,:,t)./repmat(vecnorm(E(:,:,t)),3,1);

    if mod(t,10)==0
        quiver3(X(1,:,t),X(3,:,t),X(2,:,t),E(1,:,t),E(3,:,t),E(2,:,t),'LineWidth',1,'Color',cm(t,:));
        view(2)
        daspect([1 1 1]);
        zlim([h-2,h+2]); xlim([-D-1,D+1]); ylim([-D-1,D+1]);
        drawnow();
    end


    theta = asin(E(2,1:(N-1),t));
    u = -[X(1,1:(N-1),t);zeros(1,N-1);X(3,1:(N-1),t)];  
    v= [E(1,1:(N-1),t);zeros(1,(N-1));E(3,1:(N-1),t)]; u = u./norm(u); v = v./norm(v);
    phi  = - atan2(dot([zeros(1,N-1);ones(1,N-1);zeros(1,N-1)],cross(u,v)),dot(u,v));

    ac = pi/2-pi/N; if N == 2, ac = pi/2; end
    %if max(abs(phi))> ac, disp(phi); disp(theta); disp('outwards'); return; end

    if unstable == 0 && max(abs(phi))> ac, disp('unstable azimuthal'); unstable=1; end
    if unstable == 0 && max(abs(theta))>pi/2+0.1, disp('unstable polar' ); unstable =1;  end

    if unstable ==1, disp(t); return; end
end

daspect([1 1 1]);

end

