%close all

rootdir  = 'G:\YuLin\ChR2_869661\';

savedir =  fullfile(rootdir, 'analysis_Allen\');

directionList = {'B2U', 'U2B', 'L2R', 'R2L'};

powerMaps = {};
locationMaps = {};
nanMask = {};

for i = 1:length(directionList)
        direction = directionList{i}; 
        
        fileList = dir(fullfile(savedir, sprintf('*%s*powerMap.mat', direction)));

        for ii=1:length(fileList)
            fname=strcat(savedir,fileList(ii).name);
            temp = cell2mat(struct2cell(load(fname)));
            
            %median filter
            filterRadius = 5; %note: must be odd
            temp_padded = padarray(temp, floor([filterRadius/2 filterRadius/2]), 'both'); % Pad image
            temp_paddedCol = im2col(temp_padded, [filterRadius filterRadius], 'sliding'); % Transform into columns
            temp_median = reshape(median(temp_paddedCol, 1), size(temp,1), size(temp,2)); % Find median of each column and reshape back
            
            powerMaps{i} = temp_median; 
            
            power_threshold = prctile(reshape(temp_median,[],1),90);
            powerMaps_thrMask{i} = temp_median>power_threshold;

        end
        
        fileList = dir(fullfile(savedir, sprintf('*%s*locationMap.mat', direction)));

        for ii=1:length(fileList)
            fname=strcat(savedir,fileList(ii).name);
            temp = cell2mat(struct2cell(load(fname)));
            
            %median filter
            nans = isnan(temp);
            temp_padded = padarray(temp, floor([filterRadius/2 filterRadius/2]), 'both'); % Pad image
            temp_paddedCol = im2col(temp_padded, [filterRadius filterRadius], 'sliding'); % Transform into columns
            temp_median = reshape(median(temp_paddedCol, 1), size(temp,1), size(temp,2)); % Find median of each column and reshape back
            temp_median(nans) = nan; % Set nan elements back
            
            locationMaps{i} = temp_median; 
            nanMask{i} = isnan(temp_median);
        end
end



%%
saveVariableList.azimuth_powerMap = mean(cat(3, powerMaps{3:4}), 3);
saveVariableList.altitude_powerMap = mean(cat(3, powerMaps{1:2}), 3);  


% LocationMaps currently contain: 
%   (preferred location + delay), normalised to stimulus display time for B2U/ L2R 
%   (preferred location - delay), normalised to stimulus display time for U2B/ R2L
% To obtain preferred location without delay, average B2U&U2B / L2R&R2L. 
%   Underlying math: Let x = preferred loc, d = delay and S = stimulus display,
%   (x + d)modS + (x-d)modS = (x + d + x - d)modS = (2x)modS
%   (x)modS = ( (x + d)modS + (x-d)modS )/2 
%       i.e. avg(B2U, U2B) or avg(L2R, R2L) 

saveVariableList.azimuth_locationMap = mean(cat(3, locationMaps{3:4}), 3);
saveVariableList.altitude_locationMap = mean(cat(3, locationMaps{1:2}), 3);
%%
saveVariableList.azimuth_locationMap(saveVariableList.azimuth_locationMap < 0 | saveVariableList.azimuth_locationMap > 1) = NaN;
saveVariableList.altitude_locationMap(saveVariableList.altitude_locationMap < 0 | saveVariableList.altitude_locationMap > 1) = NaN;
%%
figure; imagesc(saveVariableList.azimuth_powerMap, 'AlphaData', ~isnan(saveVariableList.azimuth_powerMap)); 
set(gca,'color',[1 1 1]);
colormap(gray); colorbar; axis square
title('Azimuth power map','Interpreter','none');

figure; imagesc(saveVariableList.altitude_powerMap, 'AlphaData', ~isnan(saveVariableList.altitude_powerMap)); 
set(gca,'color',[1 1 1]);
colormap(gray); colorbar; axis square
title('Altitude power map','Interpreter','none');

fig_alt = figure; imagesc(saveVariableList.altitude_locationMap, 'AlphaData', ~isnan(saveVariableList.altitude_locationMap)); 
set(gca,'color',[1 1 1]);
colormap(jet); colorbar; axis square
title('Altitude location map','Interpreter','none');
altitude_locationMap_fname = sprintf('%saltitude_locationFig',savedir);
print(fig_alt, altitude_locationMap_fname, '-dtiff');
    
fig_azi = figure; imagesc(saveVariableList.azimuth_locationMap, 'AlphaData', ~isnan(saveVariableList.azimuth_locationMap)); 
set(gca,'color',[1 1 1]);
colormap(jet); colorbar; axis square
title('Azimuth location map','Interpreter','none');
azimuth_locationMap_fname = sprintf('%sazimuth_locationFig',savedir);
print(fig_azi, azimuth_locationMap_fname, '-dtiff');


%%  
    
savefnameList = {'azimuth_powerMap.tif',
                 'altitude_powerMap.tif', 
                 'azimuth_locationMap.tif',
                 'altitude_locationMap.tif'};

tagstruct.ImageLength = size(saveVariableList.azimuth_powerMap,1);
tagstruct.ImageWidth = size(saveVariableList.azimuth_powerMap,2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 64;
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP; 
tagstruct.SamplesPerPixel = 1;
tagstruct.Compression = Tiff.Compression.None;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

saveVariableListCell = struct2cell(saveVariableList);
for i = 1:length(savefnameList)
        saveVar = cell2mat(saveVariableListCell(i));
        savefName = sprintf('%s%s',savedir,savefnameList{i});

        
        t = Tiff(char(savefName),'w');
        setTag(t,tagstruct);
        
        write(t, saveVar);
        close(t);
end



