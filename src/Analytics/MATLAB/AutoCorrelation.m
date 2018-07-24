%% Autocorrelation of polymer-c code

clear all;

folder      = '~/Documents/polymer-c_runs/2017_08_15_FreeSpacePolymer';
filename    = 'FreeSpacePolymerVisualize';
Nrods       = 20;

%% Read file

M = dlmread(fullfile(folder,filename));
r = zeros(size(M,1),Nrods,3);
for n = 1:Nrods
    
    r(:,n,1) = M(:,11+3*(n-1));
    r(:,n,2) = M(:,11+1+3*(n-1));
    r(:,n,3) = M(:,11+2+3*(n-1));
    
end

%% Autocorrelation
tauArray = 1:1:100;
ACF = zeros(Nrods,length(tauArray));

for tau = 1:length(tauArray)
    ACF(:,tau) = mean(mean(r(tau+1:end,:,:).*r(1:end-tau),3),1)./mean(mean(r(:,:,:).^2,3),1);
end