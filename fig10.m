function [btc] = fig10()
% Stability of circles of pushers above one or two walls
% uses the function isstable to assess the stability of a given configuration, and uses matlab function fminbnd to find the minimal |beta| that is stable 

%Parameters to be set: 

% Number of pushers (2 or 3 to have stable configurations) 
N = 2; 

% range of heights that will be explored ( L is the ratio log(epsilon)/log(h-1) ) 
L = linspace(4/3+0.1,20,15);

% if second wall, set to 2
nwalls = 1; 

% Note that the noise strength can be tuned in isstable 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Distance between neighbours
epsis = 0.2;

% returns - beta for a stable configuration and 1e3 otherwise: function to minimise 
    function [m] = stb(beta)
        if isstable(N,beta,epsis, h, nwalls)
            m = -beta;
        else
            m = 1e3;
        end
    end

btc= 0.*L;
% minimise the stable |beta |
for i = 1:length(L)
    h = 1+ exp(1/L(i)*log(epsis));
    btc(i) = fminbnd(@stb,-35,-0.6);  % range of beta explored 
end

hold on;
plot(L(btc<-0.7),btc(btc<-0.7)); 

end

