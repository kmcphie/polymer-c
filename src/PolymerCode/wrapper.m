
seriesName = 'KuhnSweepNoMem901';

argumentFile = [seriesName '.input'];

NAA = 113;
dAA = 0.3; %nm
rK = 3.4; % nm

deltaArray = (0.1:0.1:3); %nm

for iRun=1:numel(deltaArray)
    
    delta = deltaArray(iRun);
    
    N = ceil(NAA*dAA/delta);
    
    iSite = floor(0.90*N);
    
    rKRatio = rK/delta;
    
    fid = fopen(argumentFile, 'a');
    %fprintf(fid, [seriesName ' ' num2str(N) ' -1 ' num2str(rKRatio) '\n']);
    fprintf(fid, [seriesName ' ' num2str(N) ' ' num2str(iSite) ' ' num2str(rKRatio) '\n']);
    fclose(fid);
    
end
