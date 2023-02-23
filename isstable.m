function [stability] = isstable(N,beta,epsis,h,nwall)
% Returns 0 if the configuration of N squirmers of strength beta, separated
% by epsis from their neighbours and florating at a height h above a wall
% (or between two walls if nwall =2)

T = 1000;
dt = 0.5;
R = 1;
v0 =1;

if nargin== 5 && nwall == 2, wall2 = 1; elseif nwall==0, wall2 = 0; end

theta0 = 0;

E = zeros(3,N,T);

% Strenght of the noise used to assess stability

contnoise = 0.01;%0.05;
ininoise = 0.01;%0.1;
td = 1; %if 3D noise

% circle configuration
D = (epsis/2+1)./sin(pi/N);
X = [D*cos((1:1:N)*2*pi./N);h.*ones(1,N);D*sin((1:1:N)*2*pi./N)]; % position of the squirmers
E(:,:,1) = [-cos(theta0).*cos((1:1:N)*2*pi./N);sin(theta0).*ones(1,N);-cos(theta0).*sin((1:1:N)*2*pi./N)]; % normalised orientation vector

% noise in the intial orientations of the swimmers
if td == 1, E(:,:,1) = E(:,:,1)+ [2*ininoise.*(rand(1,N)-1/2);2*ininoise.*(rand(1,N)-1/2);2*ininoise.*(rand(1,N)-1/2)];
else, E(:,:,1) = E(:,:,1)+ [2*ininoise.*(rand(1,N)-1/2);0.*(rand(1,N)-1/2);2*ininoise.*(rand(1,N)-1/2)];
end

E(:,:,1)= E(:,:,1)./repmat(vecnorm(E(:,:,1)),3,1); % renormalise the orientation at each operation

for t = 2:T
    for i = 1:N
        Omi = zeros(3,1);

        % Reorientation by neighbouring swimmers
        for j = [1:i-1, i+1:N]

            epsij = norm(X(:,i)-X(:,j)) - 2;
            if epsij < 1 % the lubrication interaction occurs only below one radius separation
                eparih = (X(:,j)-X(:,i))./norm(X(:,j)-X(:,i));

                Omi = Omi + log(1./epsij)*12/(40*R)*v0*( (1+beta*dot(E(:,i,t-1), eparih) ) *cross( eparih, E(:,i,t-1)));
                Omi = Omi + log(1./epsij)*3/(40*R)*v0* ( (1+beta*dot(E(:,j,t-1), -eparih) ) *cross(-eparih, E(:,j,t-1)));

                if isnan(Omi), return; end
            end
        end

        E(:,i,t) = E(:,i,t-1)+dt.*cross(Omi,E(:,i,t-1));

    end

    % Reorientation by the wall
    if h<2 % if the height is more than 2, no lubrication interaction occurs
        for i = 1:N
            Omw = log(1./(h-1))*24/(40*R)*v0* (1 - beta * E(2,i,t-1)) *cross([0;-1;0],[E(1,i,t-1);0;E(3,i,t-1)]);
            if wall2==1
                Omw2 = log(1./(h-1))*24/(40*R)*v0* (1 + beta * E(2,i,t-1)) *cross([0;1;0],[E(1,i,t-1);0;E(3,i,t-1)]);
                E(:,i,t) = E(:,i,t) + dt.*cross(Omw2+Omw,E(:,i,t-1));
            else
                E(:,i,t) = E(:,i,t) + dt.*cross(Omw,E(:,i,t-1));
            end
        end
    end

    if td == 1
        E(:,:,t) = E(:,:,t) + contnoise.*2*(rand(3,N)-1/2).*[ones(1,N);ones(1,N);ones(1,N)];
    else
        E(:,:,t) = E(:,:,t) + contnoise.*2*(rand(3,N)-1/2).*[ones(1,N);zeros(1,N);ones(1,N)];
    end

    E(:,:,t)= E(:,:,t)./repmat(vecnorm(E(:,:,t)),3,1);

    theta = asin(E(2,:,t));
    u = -[X(1,:);zeros(1,N);X(3,:)];  v= [E(1,:,t);zeros(1,N);E(3,:,t)]; u = u./norm(u); v = v./norm(v);

    phi  = - atan2(dot([zeros(1,N);ones(1,N);zeros(1,N)],cross(u,v)),dot(u,v));

    ac = pi/2-pi/N; if N ==2, ac = pi/2; end
    if max(abs(phi))> ac, stability = 0; return; end % azimuthal instability 
    if max(abs(theta))>pi/2+0.01, stability = 0; return;  end % polar instability 
end

stability = 1;

end

