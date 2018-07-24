%% Look at polymer location distributions for Electrostatics

clear all;
close all;

%% Initialize parameters

ZETA = 1;

if(ZETA)
    combinations=64;
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
model = 2;
save = 0;
    
folder = '~/Documents/pub/lclemens/polymer-c_runs';
controlFolder = '~/Documents/polymer-c_runs';
controlSubfolder = '20180102CD3ZetaElectrostaticsDistributionControl';

saveFolder = '~/Documents/polymer-c_data_and_figures/CD3Zeta_Electrostatics/20171220';
    
switch (model)
    
    case 1
        %subfolder = '20171121CD3ZetaElectrostaticsSweep1/Dephos';
        %subfolder = '20171213CD3ZetaElectrostaticsSweep2';
        %subfolder = '20171215CD3ZetaElectrostaticsSweep3Phos1';
        %subfolder = '20171220CD3ZetaElectrostaticsSweep4Phos2';
        %subfolder = '20180110CD3ZetaElectrostaticsSweep3Phos1CorrectedPotential';
        subfolder = '20180110CD3ZetaElectrostaticsSweep4Phos2CorrectedPotential';
        filenamePrefix = 'CD3ZetaElectrostatics';
%         
%         PD = 10.^(0:0.5:2);  % parabola depth
%         PW = 10.^(-1:1:2);  % parabola width
%         WK = 10.^(-5:2:-1);  % wall parabola k
%         ER = 10.^(-2:1:2);  % repulsive energy multiplier
%         ZR = 10.^(-2:1:2);  % repulsive energy exponent multiplier


        
        PD = 10.^(0:1:2);  % parabola depth
        PW = 10.^(-1:1:2);  % parabola width
        WK = 10.^(-5:2:-1);  % wall parabola k
        ER = 10.^(0:0.5:3);  % repulsive energy multiplier
        ZR = 10.^(-1:0.5:1);  % repulsive energy exponent multiplier

    case 2
        subfolder = '20180202CD3ZetaElectrostaticsPhosphorylationCorrectedPotentialSweep7';
                

        PW = [0.1 0.3 0.5 0.7 0.9];  % parabola width
        PD = 9.*PW;  % parabola depth
        WK = 0.05;  % wall parabola k
        ER = 10.^(0:0.5:3);  % repulsive energy multiplier
        ZR = 3;  % repulsive energy exponent multiplier
        
        
        

end
    
%% Read Files

fileDoesNotExist=0;
distributionData = zeros(length(WK),length(PD),length(PW),length(ER),length(ZR),iSiteTotal+1,NBINS);

for k = 1:length(WK)
    for w = 1:length(PW)
        for d = 1:length(PD)
            for e = 1:length(ER)
                for z = 1:length(ZR)
                    
                    %% Open file and retrieve distribution data
                    switch (model)
                        case 1
                            filename = strcat(filenamePrefix,'.','PD','.',num2str(d),'PW','.',num2str(w),'WK','.',num2str(k),'ER','.',num2str(e),'ZR','.',num2str(z));      
                    end

                    if(exist(fullfile(folder,subfolder,filename))~=0)

                        M = dlmread(fullfile(folder,subfolder,filename));
                        if (ZETA)
                            ntTotal = M(1,12);
                        else
                            ntTotal = M(1,8);
                        end

                        for i=1:1:iSiteTotal+1
                            if (ZETA)
                                start = 7+9+5*iSiteTotal+3+bSiteTotal+2*iSiteTotal*bSiteTotal+iSiteTotal+(i-1)*NBINS+1;  
                            else
                                start = 3+9+5*iSiteTotal+3+bSiteTotal+2*iSiteTotal*bSiteTotal+iSiteTotal+(i-1)*NBINS+1;    
                            end
                            distributionData(k,w,d,e,z,i,:) = M(1,start:start+NBINS-1)./(ntTotal-NTCHECK);
                        end
                     else
                        fileDoesNotExist=fileDoesNotExist+1;
                        disp(filename);
                     end

                end
            end
        end
    end
end

%% Read Files
distributionDataControl = zeros(iSiteTotal+1,NBINS);

%% Open file and retrieve distribution data
controlFilename = 'DistributionControl';    

M = dlmread(fullfile(controlFolder,controlSubfolder,controlFilename));
if (ZETA)
    ntTotal = M(1,12);
else
    ntTotal = M(1,8);
end

for i=1:1:iSiteTotal+1
    if (ZETA)
        start = 7+9+5*iSiteTotal+3+bSiteTotal+2*iSiteTotal*bSiteTotal+iSiteTotal+(i-1)*NBINS+1;  
    else
        start = 3+9+5*iSiteTotal+3+bSiteTotal+2*iSiteTotal*bSiteTotal+iSiteTotal+(i-1)*NBINS+1;    
    end
    distributionDataControl(i,:) = M(1,start:start+NBINS-1)./(ntTotal-NTCHECK);
end


%% Plot histograms
lw = 1.2;
bins = 0:1:NBINS-1;
xaxis = -N+bins.*BinSize;
figureNumber = 1;
    
%for i=1:iSiteTotal+1
%for i=[1 3 5]
for i=3
    figureNumber = 100*i;
    for k=[1 2 3]
    %for k=1
        for w=[1 2 3 4]
        %for w=1
            for d=1:length(PD)
            %for d=1
                hFig = figure(figureNumber); clf; hold on;
                subfigureNumber = 1;
                for e=1:length(ER)
                %for e = 1
                    for z=1:length(ZR)
                    %for z = 1

                        plotData = reshape(distributionData(k,w,d,e,z,i,:),[NBINS,1]);
                        %subplot(length(PW),length(PD),d+length(PD)*(w-1));
                        hSub = subplot(7,5,subfigureNumber); hold on;
                        plot(xaxis,plotData,'-r','LineWidth',lw);
                        plot(xaxis,distributionDataControl(i,:),'-b','LineWidth',lw);
                        xlim([-5,10]);
                        xticks([-5 -4 -3 -2 -1 0 1 2 3 4 5]);
                        title1 = strcat('E: ', num2str(ER(e)), 'Z: ', num2str(ZR(z)));
                        %title1 = strcat('D: ', num2str(PD(d)), 'W: ', num2str(PW(w)));
                        set(gcf,'units','centimeters','position',[1,4,50,46]);
                        xlabel1 = 'z-coordinate (Kuhn lengths)';
                        ylabel1 = 'Probability';
                        %maxData(j)=max([maxData(j);plotData(:)]);
            %         
            %             if(i==iSiteTotal+1)
            %                 title1 = strcat('Distribution of Polymer Tip');
            %             else
            %                 title1 = strcat('Distribution of Tyrosine: ', num2str(i));
            %             end

                        xlabel(xlabel1,'FontName','Arial','FontSize',14);
                        ylabel(ylabel1,'FontName','Arial','FontSize',14);
                        title(title1,'FontName','Arial','FontSize',14);

                        subfigureNumber = subfigureNumber + 1;
                    end
                end
                str = {['Wall Parabola K: ',num2str(WK(k))],...
                ['Parabola Depth: ',num2str(PD(d))],...
                ['Parabola Width: ',num2str(PW(w))]};
                dim = [0.135 0.72 0.2 0.2];
                annotation(hFig,'textbox',dim,'String',str,'FitBoxToText','on');
                figureNumber = figureNumber + 1;
            end
        end
    end
end

    %% Plot energy diagrams
    if (0)
        for k=[1 3 5]
            for d=1:length(PD)
                for w=1:length(PW)
                    for i=2

                        plotData = reshape(distributionData(k,w,d,e,z,i,:),[NBINS,1]);
                        figure(((k-1)*(iSiteTotal+1)+i)*100); hold on;
                        subplot(length(PW),length(PD),d+length(PD)*(w-1));

                        plot(xaxis,KBT*log(plotData),'-r');
                        xlim([-20,20]);
                        ylabel('pNnm');
                        title(strcat('D: ', num2str(PD(d)), 'W: ', num2str(PW(w))));
                        set(gcf,'units','centimeters','position',[1,4,60,46]);

                    end
                end
            end
        end
    end

%% Save Figures

if (save)

    switch (model)
        case 1
            saveSubFolder = 'RepulsiveSweep1';
    end

    for i=15:21
        figure(i); 
        saveName = strcat('DistributioniSite',num2str(i-14));
        saveas(gcf,fullfile(saveFolder,saveSubFolder,saveName),'fig');
        saveas(gcf,fullfile(saveFolder,saveSubFolder,saveName),'epsc');
    end
    if (0)
    for i=11:17
        figure(i); 
        saveName = strcat('EnergyBarrieriSite',num2str(i-10));
        saveas(gcf,fullfile(saveFolder,saveSubFolder,saveName),'fig');
        saveas(gcf,fullfile(saveFolder,saveSubFolder,saveName),'epsc');
    end
    end

end


