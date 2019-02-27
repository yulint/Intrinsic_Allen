%close all

rootdir  = 'C:\Users\2018_Group_a\Desktop\yulin\intrinsic\GtA18\';

savedir =  fullfile(rootdir, 'analysis_Allen\');

directionList = ["B2U", "U2B", "L2R", "R2L"];

powerMaps = {};
locationMaps = {};
nanMask = {};

for i = 1:length(directionList)
        direction = directionList(i); 
        
        fileList = dir(fullfile(savedir, sprintf('*%s*powerMap.mat', direction)));

        for ii=1:length(fileList)
            fname=strcat(savedir,fileList(ii).name);
            temp = cell2mat(struct2cell(load(fname)));
            
            %median filter
            filterRadius = 1; %note: must be odd
            temp_padded = padarray(temp, floor([filterRadius/2 filterRadius/2]), 'both'); % Pad image
            temp_paddedCol = im2col(temp_padded, [filterRadius filterRadius], 'sliding'); % Transform into columns
            temp_median = reshape(median(temp_paddedCol, 1, 'omitnan'), size(temp,1), size(temp,2)); % Find median of each column and reshape back
            
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
            temp_median = reshape(median(temp_paddedCol, 1, 'omitnan'), size(temp,1), size(temp,2)); % Find median of each column and reshape back
            temp_median(nans) = nan; % Set nan elements back
            
            locationMaps{i} = temp_median; 
            nanMask{i} = isnan(temp_median);
        end
end


for i = 1:4
    
    filled_locationMaps{i} = locationMaps{i};
    
    if mod(i,2) == 1
        j = i+1;
    elseif mod(i,2) ==0
        j = i-1;
    end
    
    nonNans = ~nanMask{i} .* ~nanMask{j};
    nonNoise = powerMaps_thrMask{i} .* powerMaps_thrMask{j};
    idxs_to_regress = nonNans .* nonNoise;
    
    linear_reg_model = polyfit(squeeze(reshape(locationMaps{j}(idxs_to_regress==1),[],1)),squeeze(reshape(locationMaps{i}(idxs_to_regress==1),[],1)),1);
    predicted_i_from_j = polyval(linear_reg_model,locationMaps{j});
    
    for ii =  find(nanMask{i}==1)
        filled_locationMaps{i}(ii) = predicted_i_from_j(ii);
    
    end
    
end    

% figure; imagesc(filled_locationMaps{4}, 'AlphaData', ~isnan(filled_locationMaps{4})); 
% set(gca,'color',[1 1 1]);
% colormap(jet); colorbar; axis square
% title('nan replaced','Interpreter','none');
% figure; imagesc(filled_locationMaps{3}, 'AlphaData', ~isnan(filled_locationMaps{3})); 
% set(gca,'color',[1 1 1]);
% colormap(jet); colorbar; axis square
% title('nan replaced','Interpreter','none');
% figure; imagesc(filled_locationMaps{2}, 'AlphaData', ~isnan(filled_locationMaps{2})); 
% set(gca,'color',[1 1 1]);
% colormap(jet); colorbar; axis square
% title('nan replaced','Interpreter','none');
% figure; imagesc(filled_locationMaps{1}, 'AlphaData', ~isnan(filled_locationMaps{1})); 
% set(gca,'color',[1 1 1]);
% colormap(jet); colorbar; axis square
% title('nan replaced','Interpreter','none');


%%
saveVariableList.azimuth_powerMap = mean(cat(3, powerMaps{3:4}), 3);
saveVariableList.azimuth_locationMap = mean(cat(3, filled_locationMaps{3:4}), 3);
saveVariableList.altitude_powerMap = mean(cat(3, powerMaps{1:2}), 3);  
saveVariableList.altitude_locationMap = mean(cat(3, filled_locationMaps{1:2}), 3);



saveVariableList.azimuth_locationMap(saveVariableList.azimuth_locationMap < 0 | saveVariableList.azimuth_locationMap > 1) = NaN;
saveVariableList.altitude_locationMap(saveVariableList.altitude_locationMap < 0 | saveVariableList.altitude_locationMap > 1) = NaN;
    
    
savefnameList = [strcat(savedir, "azimuth_powerMap.tif"), 
                 strcat(savedir, "azimuth_locationMap.tif"), 
                 strcat(savedir, "altitude_powerMap.tif"), 
                 strcat(savedir, "altitude_locationMap.tif")];

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
        savefName = savefnameList(i);
        
        t = Tiff(char(savefName),'w');
        setTag(t,tagstruct);
        
        write(t, saveVar);
        close(t);
end



figure; imagesc(saveVariableList.azimuth_locationMap, 'AlphaData', ~isnan(saveVariableList.azimuth_locationMap)); 
set(gca,'color',[1 1 1]);
colormap(jet); colorbar; axis square
title('Azimuth location map','Interpreter','none');

figure; imagesc(saveVariableList.azimuth_powerMap, 'AlphaData', ~isnan(saveVariableList.azimuth_powerMap)); 
set(gca,'color',[1 1 1]);
colormap(gray); colorbar; axis square
title('Azimuth power map','Interpreter','none');

figure; imagesc(saveVariableList.altitude_locationMap, 'AlphaData', ~isnan(saveVariableList.altitude_locationMap)); 
set(gca,'color',[1 1 1]);
colormap(jet); colorbar; axis square
title('Altitude location map','Interpreter','none');

figure; imagesc(saveVariableList.altitude_powerMap, 'AlphaData', ~isnan(saveVariableList.altitude_powerMap)); 
set(gca,'color',[1 1 1]);
colormap(gray); colorbar; axis square
title('Altitude power map','Interpreter','none');
