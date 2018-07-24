%% Function to read in files

% WOULD all of this be easier with a structure?

%function [] = Funct_EffConcReadFiles(savename,Nrods,irLigand,folder,subfolderprefix,dataThin)
function [] = Funct_EffConcReadFiles(savename,dataFolder,Nrods,irLigand,dataThin)

% Initialize data cell arrays
% size: length of Nrods, length of irLigand, membrane on/off
ligandBoundrData        = cell(size(Nrods,2),size(irLigand,2),2);
ligandOffrData          = cell(size(Nrods,2),size(irLigand,2),2);
ligandBoundcenterData   = cell(size(Nrods,2),size(irLigand,2),2);
ligandOffcenterData     = cell(size(Nrods,2),size(irLigand,2),2);
occlusionData           = cell(size(Nrods,2),size(irLigand,2),2);

    for membrane=0:1:1

        % construct appropriate folder/filenames
        % REALLY need to have this completely separate - shouldn't need to
        % modify this
        if (membrane)
            subfolderprefix = 'HalfSpace_';
        else
            subfolderprefix = 'FreeSpace_';
        end
% 
%         foldername1 = strcat(folder,'June072017SigmaRSweep',subfoldersuffix,'LigandOn');
%         foldername2 = strcat(folder,'June072017SigmaRSweep',subfoldersuffix,'LigandOff');
% 
%         subfolder1  = strcat(subfolderprefix,subfoldersuffix, 'LigandOn');
%         subfolder2  = strcat(subfolderprefix,subfoldersuffix, 'LigandOff');

          foldername1 = dataFolder;
          subfolder1  = strcat(subfolderprefix,'LigandBound');
          
          foldername2 = dataFolder;
          subfolder2  = strcat(subfolderprefix,'LigandOff');


        for n = 1:length(Nrods)
            
            disp(Nrods(n));
            
            for irL = 1:length(irLigand)
                
                disp(irLigand(irL));
                % clear relevant variables between runs
                clear M1 M2 r1 r2 ligandCenter1 ligandCenter2 occlusion;

                % create filenames
                filename1 = strcat(subfolder1,'.N.',num2str(Nrods(n)),'.R.',num2str(irLigand(irL)));
                filename2 = strcat(subfolder2,'.N.',num2str(Nrods(n)),'.R.',num2str(irLigand(irL)));

                % read in data
                M1 = dlmread(fullfile(foldername1,subfolder1,filename1));
                M2 = dlmread(fullfile(foldername2,subfolder2,filename2));

                % match file sizes, either by trimming or expanding dataset
                if (size(M1,1) > size(M2,1))
                    if (dataThin)
                        % if trimming data down to size, randomly select rows
                        % to remove from larger dataset
                        MScrap  = M1;
                        cutRows = randi([1 size(M1,1)],size(M2,1),1);
                        M1      = MScrap(cutRows(:),:);
                    else
                        % if expanding data, randomly add data back into
                        % dataset, then double both datasets by randomly
                        % choosing and attaching more data from each dataset
                        MScrap  = M2;
                        addRows = randi([1 size(M2,1)],abs(size(M1,1)-size(M2,1)),1);
                        M2      = [M2; MScrap(addRows(:),:);];

                        % double the size of the data
                        addRows1 = randi([1 size(M1,1)],size(M1,1),1);
                        addRows2 = randi([1 size(M2,1)],size(M2,1),1);
                        M1       = [M1; M1(addRows1(:),:)]; % removed second semi colon, not sure if it did anything
                        M2       = [M2; M2(addRows2(:),:)];
                    end
                elseif (size(M1,1) < size(M2,1))
                    if (dataThin)
                        MScrap  = M2;
                        cutRows = randi([1 size(M2,1)],size(M1,1),1);
                        M2      = MScrap(cutRows(:),:);
                    else
                        MScrap  = M1;
                        addRows = randi([1 size(M1,1)],abs(size(M2,1)-size(M1,1)),1);
                        M1      = [M1; MScrap(addRows(:),:)];

                        % double the size of the data
                        addRows1 = randi([1 size(M1,1)],size(M1,1),1);
                        addRows2 = randi([1 size(M2,1)],size(M2,1),1);
                        M1       = [M1; M1(addRows1(:),:)];
                        M2       = [M2; M2(addRows2(:),:)];
                    end

                end

                % collect data for end point location and ligand center location
                for i=1:3
                    r1(:,i)             = M1(:, 11+(i-1));
                    r2(:,i)             = M2(:, 11+(i-1));
%                     ligandCenter1(:,i)  = M1(:, 7+3+(i-1));
                    ligandCenter1(:,i)  = M1(:, 11+3+(i-1));
                    ligandCenter2(:,i)  = M2(:, 11+3+(i-1));
                end

                % collect occlusion data from second data set
                occlusion(:) = M2(:,11+3+3+1);

                % create cell array of data
                ligandBoundrData{n,irL,membrane+1}      = r1;
                ligandOffrData{n,irL,membrane+1}        = r2;
                ligandBoundcenterData{n,irL,membrane+1} = ligandCenter1;
                ligandOffcenterData{n,irL,membrane+1}   = ligandCenter2;
                occlusionData{n,irL,membrane+1}         = occlusion;

            end
        end
    end
    
disp(savename);   
if(exist(savename)==2)
    disp('Exists!');
else
    save(savename,'ligandBoundrData','ligandOffrData','ligandBoundcenterData','ligandOffcenterData','occlusionData');
end


end