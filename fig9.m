function [] = fig9
% Stability of circles of pullers above one or two walls
% uses the function isstable to assess the stability of a given configuration, and uses matlab function fminbnd to find the minimal |beta| that is stable 

%Parameters to be set: 

% Range for the number of pullers (set above 5 to have stable configurations) 
NN = 5:6;%15; 

% range of heights that will be explored ( L is the ratio log(epsilon)/log(h-1) )  
LL =  [.3 .5 1 1.5]; % linspace(4/3+0.1,20,15);

% if second wall, set to 2 (to get fig11) 
nwalls = 1; 

% Note that the noise strength can be tuned in isstable 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% returns beta for a stable configuration and 1e3 otherwise: function to minimise 
    function [m] = stb(beta)
        if isstable(N,beta,epsis,h,nwalls)
            m = beta;
        else
            m = 1e3;
        end
    end

epsis = 0.5;

btc= 0.*NN;

cm = colormap(hot(length(LL)+4));

% minimise the stable beta
for j = 1:length(LL)
    h = 1 + exp(1/LL(j)*log(epsis));
     % j % use as countdown 
     for i = 1:length(NN)
         N=NN(i);
         btc(i) = fminbnd(@stb,0.01,18); % range of explored beta 
     end 
    hold on;
    plot(NN, 5 .*LL(j)*sin(pi./NN)./(8-5.*LL(j).*sin(pi./NN).^2),':','LineWidth',1.2,'Color',cm(end-j-3,:));
    plot(NN,btc,'.','MarkerSize',16,'Color',cm(end-j-3,:));
end
end

