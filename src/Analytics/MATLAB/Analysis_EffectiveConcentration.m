%% Analysis_EffectiveConcKernel
% Sigma(r) computations - Edited to remove parameter dependence
% Cell array version

clear all;
close all;

%% Initialize parameters
NRODS       = [1];                  % number of segments
delta       = 1;                    % Kuhn length
rho         = 0:1:10;               % separation distance between two polymers
irLigand    = [1];                  % ligand size
cutoff      = 2;                    % cutoff for ligand being 'close' to binding site, in Kuhn lengths
%stDevDex   = (5.6/0.3)/(sqrt(3));  % divide by 0.3 to put in kuhn lengths, divide by sqrt(3) to get each individual dimension stdev
                                        %set to zero to remove influence of dextran
stDevDex    = 0;
dataThin    = 0;
NTCHECK     = 200000;
saveTF      = 0;

%% Folder and file names
folder          = '~/Documents/polymer-c_runs/2017_08_15_EffectiveConcentrationSingleSegment';
%subfolderprefix = 'SPROcclusion';

% filename with appropriate data, or place to save appropriate data
dataFilename    = strcat('2017_08_22_SPRData','_N_',num2str(NRODS(1)),'-',num2str(NRODS(end)),...
                         '_irL_',num2str(irLigand(1)),'-',num2str(irLigand(end)));

savefolder      = '~/Google Drive/polymer-c/polymer-c_data_and_figures/SurfaceEffects/EffectiveConcentrationKernel/UpdatedJune072017';
 
dataFolder = '~/Documents/polymer-c_runs/2017_08_15_EffectiveConcentrationSingleSegment';
%% Read in data

% PLAY WITH THIS - include full file path, possibly 'file' restriction
if( exist(strcat(dataFilename,'.mat'),'file') == 2 )
    load(dataFilename);
else
%     Funct_EffConcReadFiles(dataFilename,NRODS,irLigand,folder,subfolderprefix,dataThin);
    Funct_EffConcReadFiles(dataFilename,dataFolder,NRODS,irLigand,dataThin);
    load(dataFilename);
end


%% Initialize storage arrays
            
totalClose          = zeros(length(NRODS),length(irLigand),2,length(rho));
totalAbleToBind     = zeros(length(NRODS),length(irLigand),2,length(rho));
fraction            = zeros(length(NRODS),length(irLigand),2,length(rho));
sigma               = zeros(length(NRODS),length(irLigand),2,length(rho));
sigmaError          = zeros(length(NRODS),length(irLigand),2,length(rho));
sigmaNondim         = zeros(length(NRODS),length(irLigand),2,length(rho));
sigmaErrorNondim    = zeros(length(NRODS),length(irLigand),2,length(rho));
cN                  = zeros(length(NRODS),length(irLigand),2);
sigmaTheory         = zeros(length(NRODS),200);

%% Calculate effective concentration    
for n = 1:length(NRODS)
    for membrane=0:1:1
        for irL = 1:length(irLigand)
            clear ligandCenter1 occlusion availability;
            
            % set up vector of bound ligand locations
            ligandCenter1 = ligandBoundcenterData{n,irL,membrane+1};
            occlusion     = occlusionData{n,irL,membrane+1};
            availability  = ~occlusion; % need availability, not occlusion

                for r = 1:length(rho)
                    clear ligandCenter2shifted ligandOffrDatashifted distance distanceClose bindingAvailability;
                    
                    %%%%%% Calculate shifted polymer data
                    
                    % Initialize shifted data as original data

                    ligandCenter2shifted  = ligandOffcenterData{n,irL,membrane+1};
                    ligandOffrDatashifted = ligandOffrData{n,irL,membrane+1};
                    
                    % shift unbound polymer data by rho(r)in the x
                    % direction

                    ligandCenter2shifted(:,1)  = ligandCenter2shifted(:,1)+rho(r);
                    ligandOffrDatashifted(:,1) = ligandOffrDatashifted(:,1)+rho(r);

                    % create random matrix of Dextran shifts - only use if in
                    % free-space

                    if (~membrane)
                        
                        dextranShift1        = normrnd(0,stDevDex,[size(ligandCenter1,1),size(ligandCenter1,2)]);
                        ligandCenter1        = ligandCenter1+dextranShift1;
                        
                        dextranShift2        = normrnd(0,stDevDex,[size(ligandCenter2shifted,1),size(ligandCenter2shifted,2)]);
                        ligandCenter2shifted = ligandCenter2shifted+dextranShift2;
                        % SHOULD also shift r data in case want it
                    end
                    
                    
                    
                    %% CONSIDER INSERTING CALL TO VISUALIZER HERE
                    
                    %%
                    %%%%%% Calculate effective concentration:
                    
                    % calculate distance between bound ligand center and
                    % binding site center

                    distance = sqrt((ligandCenter1(:,1)-ligandCenter2shifted(:,1)).^2 +...
                                    (ligandCenter1(:,2)-ligandCenter2shifted(:,2)).^2 +...
                                    (ligandCenter1(:,3)-ligandCenter2shifted(:,3)).^2);
                    
                    % determine when ligand is close enough to binding site
                    % based on cutoff distance

                    distanceClose = (distance<cutoff);

                    % calculate how many instances are close enough
                    totalClose(n, irL, membrane+1, r) = size(find(distanceClose==1),1);
                    
                    % determine which instances are both close enough AND
                    % the binding site is available
                    bindingAvailability = distanceClose.*availability';
                    
                    % find total both close and unoccluded
                    totalAbleToBind(n, irL, membrane+1, r) = size(find(bindingAvailability),1);
                    disp(totalAbleToBind(1,1,1,:));
                    % find fraction of instances close and unoccluded
                    % compared to total instances
                    fraction(n, irL, membrane+1, r) = size(find(bindingAvailability),1)/...
                                                          size(bindingAvailability,1);

                    % calculate error for fraction
                    BernoulliError = sqrt( fraction(n, irL, membrane+1, r)*(1-fraction(n, irL, membrane+1, r) )/...
                                     size(bindingAvailability,1) );
                    
                    % find effective concentration
                        % divide fraction (close and unoccluded) by volume
                        % (close enough cutoff sphere)
                    sigma(n, irL, membrane+1, r) = fraction(n, irL, membrane+1, r)/(4/3*pi*(cutoff)^3);

                    % calculate error for effective concentration
                    sigmaError(n, irL, membrane+1, r) = BernoulliError / (4/3*pi*(cutoff)^3);
                end
        end
    end
    

    %% Calculate magnitude of concentration ehancement from sigma at r = 0
    % where does the 2*pi come from?  I expect it should be 4*pi*N?
    %cN  = (4*pi/3)^(3/2)*(1/sigma(1));
    
    % cN(NRODS, IRLIGAND, MEMBRANE)
    cN(n, :, :) = (4*pi*delta^2*NRODS(n)/3)^(3/2)*sigma(n, :, :, 1);
    
    %% Sigma(r) nondimensionalized

    sigmaNondim(n,:,:,:)      = sigma(n,:,:,:).*delta^3.*NRODS(n)^(3/2);
    sigmaErrorNondim(n,:,:,:) = sigmaError(n,:,:,:).*delta^3.*NRODS(n)^(3/2);

    %% Calculate sigmaTheory
    % Theoretical curve for Worm-like Chain model
    % COMPARE to equation in Goyette et al. 2017 (probably not the same....)
    rhoArray             = linspace(0,max(rho),200);
    sigmaTheory(n,:)     = ( 3/(4*pi*delta^2*NRODS(n)) )^(3/2) * exp( (-3*rhoArray.^2)/(4*NRODS(n)*delta^2) );

end

%% Calculate ratio of sigma half space to free space, with error

ratioHalftoFreeSpace = sigma(:,:,2,:)./sigma(:,:,1,:);
ratioError           = abs(sigmaNondim(:,:,2,:)./sigmaNondim(:,:,1,:)).*...
                        sqrt((sigmaErrorNondim(:,:,1,:)./sigmaNondim(:,:,1,:)).^2+...
                             (sigmaErrorNondim(:,:,2,:)./sigmaNondim(:,:,2,:)).^2);

%% Save variables for plotting

saveFilename = strcat(dataFilename,'_EffConcVars');
save(saveFilename,'NRODS','irLigand','delta','rho','sigma','sigmaError','rhoArray','sigmaTheory','sigmaNondim','sigmaErrorNondim','ratioHalftoFreeSpace','ratioError');
    
%% Plot effective concentration

Funct_EffConcPlots(saveFilename,savefolder,saveTF);
    

