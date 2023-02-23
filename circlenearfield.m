function [] = circlenearfield()
% This function explores the dynamics of a circle of squirmers in the lubrication regime 
% It provides a demonstration of how the function isstable computed the stability of the cluster 
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
td = 1;                 %if 3D noise (as opposed to in-plane for phi only) 

X = zeros(3,N,T);
E = zeros(3,N,T);

cm = colormap(hot(floor(T*1.2)));


% circle configuration
D = (epsis/2+1)./sin(pi/N);
X(:,:,1) = [D*cos((1:1:N)*2*pi./N);h.*ones(1,N);D*sin((1:1:N)*2*pi./N)];
E(:,:,1) = [-cos(theta0).*cos((1:1:N)*2*pi./N);sin(theta0).*ones(1,N);-cos(theta0).*sin((1:1:N)*2*pi./N)];


% noise
if td == 1, E(:,:,1) = E(:,:,1)+ [2*ininoise.*(rand(1,N)-1/2);2*ininoise.*(rand(1,N)-1/2);2*ininoise.*(rand(1,N)-1/2)];
else %E(:,:,1) = E(:,:,1) + [zeros(1,N); 0.1,-0.1,0.1,-0.1,0.1, -0.1,0.1,-0.1,0.1,-0.1; zeros(1,N)];
    E(:,:,1) = E(:,:,1) + [2*ininoise.*(rand(1,N)-1/2);0.*(rand(1,N)-1/2);2*ininoise.*(rand(1,N)-1/2)];
end

E(:,:,1)= E(:,:,1)./repmat(vecnorm(E(:,:,1)),3,1);

figure
clf
hold on;

quiver3(X(1,:,1),X(3,:,1),X(2,:,1),E(1,:,1),E(3,:,1),E(2,:,1),'LineWidth',1,'Color',cm(1,:));
drawnow();
xlabel('x'); ylabel('z');


unstable =0;

for t = 2:T
    for i = 1:N
        Omi = zeros(3,1);
        for j = [1:i-1, i+1:N]

            epsij = norm(X(:,i,t-1)-X(:,j,t-1)) - 2;
            if epsij < 1

                %disp([i,j])
                eparih = (X(:,j,t-1)-X(:,i,t-1))./norm(X(:,j,t-1)-X(:,i,t-1));

                Omi = Omi + log(1./epsij)*12/(40*R)*v0*( (1+beta*dot(E(:,i,t-1), eparih) ) *cross( eparih, E(:,i,t-1)));
                Omi = Omi + log(1./epsij)*3/(40*R)*v0* ( (1+beta*dot(E(:,j,t-1), -eparih) ) *cross(-eparih, E(:,j,t-1)));

                %                 if i == N
                %                     j
                %                     dot(E(:,i,t-1), -eparih)
                %                     cos(asin(E(2,i,t-1)))*sin(pi/N)
                %                     cross(eparih, E(:,i,t-1))
                %                     cross(-eparih, E(:,j,t-1))
                %                     sin(asin(E(2,i,t-1)))*sin(pi/N)
                %
                %                     if j ==8, return; end
                %                 end
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

    if td == 1
        E(:,:,t) = E(:,:,t) + contnoise.*2*(rand(3,N)-1/2).*[ones(1,N);ones(1,N);ones(1,N)];
    else
        E(:,:,t) = E(:,:,t) + contnoise.*2*(rand(3,N)-1/2).*[ones(1,N);zeros(1,N);ones(1,N)];
    end

    E(:,:,t)= E(:,:,t)./repmat(vecnorm(E(:,:,t)),3,1);

    if mod(t,10)==0
        quiver3(X(1,:,t),X(3,:,t),X(2,:,t),E(1,:,t),E(3,:,t),E(2,:,t),'LineWidth',1,'Color',cm(t,:));
        view(2)
        daspect([1 1 1]);
        zlim([h-2,h+2]); xlim([-D-1,D+1]); ylim([-D-1,D+1]);
        drawnow();
    end

    theta = asin(E(2,:,t));
    u = -[X(1,:,t);zeros(1,N);X(3,:,t)];  v= [E(1,:,t);zeros(1,N);E(3,:,t)]; u = u./norm(u); v = v./norm(v);
    phi  = - atan2(dot([zeros(1,N);ones(1,N);zeros(1,N)],cross(u,v)),dot(u,v));

    ac = pi/2-pi/N; if N ==2, ac = pi/2; end
    % if max(abs(phi))> ac, disp(phi); disp(theta); disp('outwards'); return; end

    if unstable == 0 && max(abs(phi))> ac, disp('unstable azimuthal'); unstable=1; end
    if unstable == 0 && max(abs(theta))>pi/2+0.1, disp('unstable polar' ); unstable =1;  end

    if unstable ==1, disp(t); return; end
end

disp(['\theta _{max} = ', num2str(mean(abs(theta)))]);
disp(['\delta \phi _{max} = ',num2str(max(abs(phi)))]);

daspect([1 1 1]);


end

