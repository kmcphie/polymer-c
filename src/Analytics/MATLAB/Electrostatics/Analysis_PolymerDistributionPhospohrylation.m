%% Look at polymer location distributions for Electrostatics

clear all;
close all;

%% Initialize parameters


ZETA = 1;

if(ZETA)
    combinations=1;
    N=113;
    iSiteTotal = 6;
else
    combinations = 4;
    N = 57;
    iSiteTotal = 2;
end

NTCHECK = 200000;
NBINS = 3000;
KBT = 4.14;
BinSize = 2*N/NBINS;
bSiteTotal = 0;
model = 10;
save = 0;

    
%folder = '~/Documents/pub/lclemens/polymer-c_runs';
folder = '~/Documents/polymer-c_runs';
    
switch (model)
    
    case 1
        subfolder = 'Apr262017CD3ZetaPhosphorylationHardwallPiecewiseBasicsY';
        filenamePrefix = 'CD3ZetaPhosphorylationHardwallPiecewiseBasicsY';
        depthName = 'ParabolaDepth';
        widthName = 'ParabolaWidth';
        
        depth = 1; %ParabolaDepth
        width = 1; %ParabolaCenter  
        wallParabolaK =1;
        
    case 2
        subfolder = 'Apr262017CD3ZetaPhosphorylationHardwallPiecewiseYOnly';
        filenamePrefix = 'CD3ZetaPhosphorylationHardwallPiecewiseYOnly';
        depthName = 'ParabolaDepth';
        widthName = 'ParabolaWidth';
        
        depth = 1; %ParabolaDepth
        width = 1; %ParabolaCenter  
        wallParabolaK =1;
        
    case 3
        subfolder = 'Apr262017CD3EpsilonPhosphorylationHardwallPiecewiseBasicsYDepthSweepFine';
        filenamePrefix = 'CD3EpsilonPhosphorylationHardwallPiecewiseBasicsY';
        depthName = 'ParabolaDepth';
        widthName = 'ParabolaWidth';
        
        %depth = 1; %ParabolaDepth
        depth = [4 6 8 10];
        width = 1; %ParabolaCenter  
        wallParabolaK =1;
        
    case 4
        subfolder = 'Apr262017CD3EpsilonPhosphorylationHardwallPiecewiseYOnlyDepthSweepFine';
        filenamePrefix = 'CD3EpsilonPhosphorylationHardwallPiecewiseYOnly';
        depthName = 'ParabolaDepth';
        widthName = 'ParabolaWidth';
        
        %depth = 1; %ParabolaDepth
        depth = [4 6 8 10];
        width = 1; %ParabolaCenter  
        wallParabolaK =1;
        
    case 5
        subfolder = 'Apr262017CD3ZetaPhosphorylationSoftwallPiecewiseBasicsY';
        filenamePrefix = 'CD3ZetaPhosphorylationSoftwallPiecewiseBasicsY';
         subfolderControl = 'May232017CD3ZetaSoftwallPiecewiseBasicsYPhosphorylationControl';
        filenamePrefixControl = 'CD3ZetaSoftwallPiecewiseBasicsYPhosphorylationControl';
        depthName = 'ParabolaDepth';
        widthName = 'ParabolaWidth';
        wallParabolaKName = 'WallParabolaK';
        
        depth = 1; %ParabolaDepth
        width = 1; %ParabolaCenter  
        wallParabolaK =1;
        
    case 6
        subfolder = 'Apr262017CD3ZetaPhosphorylationSoftwallPiecewiseYOnly';
        filenamePrefix = 'CD3ZetaPhosphorylationSoftwallPiecewiseYOnly';
        subfolderControl = 'May232017CD3ZetaSoftwallPiecewiseBasicsYPhosphorylationControl';
        filenamePrefixControl = 'CD3ZetaSoftwallPiecewiseBasicsYPhosphorylationControl';
        depthName = 'ParabolaDepth';
        widthName = 'ParabolaWidth';
        wallParabolaKName = 'WallParabolaK';
        
        depth = 1; %ParabolaDepth
        width = 1; %ParabolaCenter  
        wallParabolaK =1;

    case 7
        subfolder = 'Apr262017CD3EpsilonPhosphorylationSoftwallPiecewiseBasicsYDepthWallWidthSweepFine';
        filenamePrefix = 'CD3EpsilonPhosphorylationSoftwallPiecewiseBasicsY';
        depthName = 'ParabolaDepth';
        widthName = 'ParabolaWidth';
        wallParabolaKName = 'WallParabolaK';
        
        %depth = 1; %ParabolaDepth
        depth = [4 6 8 10];
        width = 1; %ParabolaCenter  
        wallParabolaK =[1 3.1 10 31];
        
    case 8
        subfolder = 'Apr262017CD3EpsilonPhosphorylationSoftwallPiecewiseYOnlyDepthWallWidthSweepFine';
        filenamePrefix = 'CD3EpsilonPhosphorylationSoftwallPiecewiseYOnly';
        depthName = 'ParabolaDepth';
        widthName = 'ParabolaWidth';
        wallParabolaKName = 'WallParabolaK';
        
        %depth = 1; %ParabolaDepth
        depth = [4 6 8 10];
        width = 1; %ParabolaCenter  
        wallParabolaK =[1 3.1 10 31];
        
     case 9
        subfolder = 'May232017CD3ZetaSoftwallPiecewiseBasicsYPhosphorylation';
        filenamePrefix = 'CD3ZetaSoftwallPiecewiseBasicsYPhosphorylation';
        subfolderControl = 'May232017CD3ZetaSoftwallPiecewiseBasicsYPhosphorylationControl';
        filenamePrefixControl = 'CD3ZetaSoftwallPiecewiseBasicsYPhosphorylationControl';
        depthName = 'ParabolaDepth';
        widthName = 'ParabolaWidth';
        wallParabolaKName = 'WallParabolaK';
        
        %depth = 1; %ParabolaDepth
        depth = 10;
        width = 1; %ParabolaCenter  
        wallParabolaK = 0.01;
        
      case 10
        subfolder = 'July022017CD3ZetaSoftwallPiecewiseBasicsOnly';
        filenamePrefix = 'CD3ZetaSoftwallPiecewiseBasicsOnly';
        subfolderControl = 'July022017CD3ZetaSoftwallPiecewiseBasicsOnlyMutated';
        filenamePrefixControl = 'CD3ZetaSoftwallPiecewiseBasicsOnly';
        depthName = 'ParabolaDepth';
        widthName = 'ParabolaWidth';
        wallParabolaKName = 'WallParabolaK';
        
        %depth = 1; %ParabolaDepth
        depth = 31.6;
        width = 1; %ParabolaCenter  
        wallParabolaK = 10;

end
    

for k=1:length(wallParabolaK)
for d=1:length(depth)
%% Read Files



distributionData = zeros(iSiteTotal+1,combinations,NBINS);


%filename = strcat(filenamePrefix,'ParabolaDepth.',num2str(d),'.ParabolaWidth.1.WallParabolaK.',num2str(k),'.cat.txt');
%filename = strcat(filenamePrefix,'ParabolaDepth.',num2str(d),'.cat.txt');
%filename = strcat(filenamePrefix,'.cat.txt');
filename = 'CD3ZetaPhosphorylationBasicsOnly';
%filenameControl = strcat(filenamePrefixControl,'ParabolaDepth.1.ParabolaWidth.1.WallParabolaK.1');
filenameControl = 'CD3ZetaPhosphorylationBasicsOnly';
       M = dlmread(fullfile(folder,subfolder,filename));
       MControl = dlmread(fullfile(folder,subfolderControl,filenameControl));
            
       

       
       for j=1:combinations
            if (ZETA)
                ntTotal = M(j,12);
            else
                ntTotal = M(j,8);
            end
        
            for i=1:1:iSiteTotal+1
                
                        if (ZETA)
                            start = 7+9+4*iSiteTotal+3+bSiteTotal+(i-1)*NBINS+1;  
                        else
                            start = 3+9+4*iSiteTotal+3+bSiteTotal+(i-1)*NBINS+1;    
                        end

                distributionData(i,j,:) = M(j,start:start+NBINS-1)./(ntTotal-NTCHECK);
                
            end
       end
       
       for i=1:1:iSiteTotal+1
           
                   if (ZETA)
                        ntTotal = MControl(1,12);
                        start = 7+9+4*iSiteTotal+3+bSiteTotal+(i-1)*NBINS+1;  
                    else
                        start = 3+9+4*iSiteTotal+3+bSiteTotal+(i-1)*NBINS+1;  
                        ntTotal = MControl(1,8);
                    end
                distributionDataControl(i,:) = MControl(1,start:start+NBINS-1)./(ntTotal-NTCHECK);
                
       end

        




%% Plot histograms

bins = 0:1:NBINS-1;
xaxis = -N+bins.*BinSize;

colors = parula(20);
%seq=[1 7 22 42 57 63 64];
%seq=[1 2 3 4];
seq = 1;
for j=1:length(seq)
    
    maxData(j)=0;
    for i=1:iSiteTotal+1
        plotData = reshape(distributionData(i,seq(j),:),[NBINS,1]);
        figure(k); hold on;
        subplot(2,4,i); 
        %subplot(length(depth),iSiteTotal+1,(iSiteTotal+1)*(d-1)+i);
        plot(xaxis,plotData,'Color',colors(j*2,:),'LineWidth',1.2);
        xlim([-5,5]);
        xticks(-5:1:5);
        %set(gcf,'units','centimeters','position',[1,4,50,20]);
        xlabel1 = 'z-coordinate (Kuhn lengths)';
        ylabel1 = 'Probability';
        maxData(j)=max([maxData(j);plotData(:)]);
        
        if(i==iSiteTotal+1)
            title1 = 'Distribution of Polymer Tip';
        elseif (i==iSiteTotal+2)
            title1 = 'Distribution of Non-basic, non-tyrosine';
        else
            title1 = strcat('  Distribution of Tyrosine: ', num2str(i));
        end
 
        xlabel(xlabel1,'FontName','Arial','FontSize',18);
        ylabel(ylabel1,'FontName','Arial','FontSize',18);
        title(title1,'FontName','Arial','FontSize',18);
        
        
        set(gcf,'units','centimeters','position',[1,4,60,46]);
    end
    
    

end


for i=1:iSiteTotal+1
    subplot(2,4,i);
    %subplot(length(depth),3,(iSiteTotal+1)*(d-1)+i);
    %plot(xaxis,distributionDataControl(i,:),'Color',colors((j+1)*2,:),'LineWidth',1.2);
    line([0 0],[0 1],'LineWidth',1.0,'LineStyle','--','Color','k');
    %legend('None', 'Y2','Y1','Both','Membrane Location');
    legend('None', '100000','110000','111000','111100','111110','111111','No BRs or Y');
    ylim([0 max(maxData)]);
end



figure(k);
suptitle('Polymer distribution as phosphorylation state changes');



%%
if(1)
%seq = [1 2 8 23 43 58 64];
%seq = [1 7 21 37 49 61 64];
seq = 1;
for j=1:length(seq)
    
    maxData(j)=0;
    for i=1:iSiteTotal+1
        plotData = reshape(distributionData(i,seq(j),:),[NBINS,1]);
        figure(i+1); hold on;
        %subplot(2,4,i); 
        plot(xaxis,plotData,'Color',colors(j*2,:),'LineWidth',1.2);
        xlim([-10,10]);
        xticks(2*[-5 -4 -3 -2 -1 0 1 2 3 4 5]);
        %set(gcf,'units','centimeters','position',[1,4,50,20]);
        xlabel1 = 'z-coordinate (Kuhn lengths)';
        ylabel1 = 'Probability';
        maxData(j)=max([maxData(j);plotData(:)]);
        
        if(i==iSiteTotal+1)
            title1 = 'Distribution of Polymer Tip';
        elseif (i==iSiteTotal+2)
            title1 = 'Distribution of Non-basic, non-tyrosine';
        else
            title1 = strcat('  Distribution of Tyrosine: ', num2str(i));
        end
 
        xlabel(xlabel1,'FontName','Arial','FontSize',18);
        ylabel(ylabel1,'FontName','Arial','FontSize',18);
        title(title1,'FontName','Arial','FontSize',18);
        
        
        set(gcf,'units','centimeters','position',[1,1,40,35]);
    end
    %color=color+5;
    
    

end

for i=1:iSiteTotal+1
    figure(i+1);
    %subplot(2,4,i);
    plot(xaxis,distributionDataControl(i,:),'Color',colors((7)*2+1,:),'LineWidth',1.2);
    line([0 0],[0 1],'LineWidth',1.0,'LineStyle','--','Color','k');
    %legend('None', '100000','101000','101010','101011','111011','111111','No BRs or Y');
    legend('WT Basic Residues','Mutated Basic Residues');
    ylim([0 max(maxData)]);
end



%figure(2);
%suptitle('Polymer distribution as phosphorylation state changes');
    
end

%% Plot energy diagrams
if (0)
    for j=[1 64]
            for i=1:iSiteTotal+1

                plotData = reshape(distributionData(i,j,:),[NBINS,1]);
                figure(2); hold on;
                subplot(3,3,i);

                plot(xaxis,KBT*log(plotData),'-r');
                xlim([-20,20]);
                ylabel('pNnm');
                title(strcat('iSite: ',num2str(i)));
                set(gcf,'units','centimeters','position',[1,4,60,46]);

            end
    end
end
end

%% Save Figures

if (save)
    
    if(model== 1 | model==2 | model==5 | model==6 | model==9 )
        saveFolder = '~/Documents/polymer-c_data_and_figures/CD3Zeta_ElectrostaticsPhosphorylation';
    else
        saveFolder = '~/Documents/polymer-c_data_and_figures/CD3Epsilon_ElectrostaticsPhosphorylation';
    end
    
    saveSubFolder = filenamePrefix;
    
    for i=1
        figure(k); 
        %saveName = 'iSiteDistributionDepth31Width1';
         saveName = 'iSiteDistribution123456';
        saveas(gcf,fullfile(saveFolder,saveSubFolder,saveName),'fig');
        saveas(gcf,fullfile(saveFolder,saveSubFolder,saveName),'epsc');
    end
    
    if (1)
    for i=2
        figure(i); 
        saveName = 'iSiteDistribution135624';
        saveas(gcf,fullfile(saveFolder,saveSubFolder,saveName),'fig');
        saveas(gcf,fullfile(saveFolder,saveSubFolder,saveName),'epsc');
    end
    end
    
    if(0)
    for i=11:17
        figure(i); 
        saveName = strcat('EnergyBarrieriSite',num2str(i-10));
        saveas(gcf,fullfile(saveFolder,saveSubFolder,saveName),'fig');
        saveas(gcf,fullfile(saveFolder,saveSubFolder,saveName),'epsc');
    end
    end
    
end

end


