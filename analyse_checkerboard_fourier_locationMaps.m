rootdir  = 'C:\Users\2018_Group_a\Desktop\yulin\intrinsic\';

savedir =  fullfile(rootdir, 'analysis_Allen\');

directionList = ["B2U", "U2B", "L2R", "R2L"];

powerMaps = {};
locationMaps = {};

for i = 1:length(directionList)
        direction = directionList(i);
        
        fileList = dir(fullfile(savedir, sprintf('*%s*powerMap.mat', direction)));

        for ii=1:length(fileList)
            fname=strcat(savedir,fileList(ii).name);
            temp = cell2mat(struct2cell(load(fname)));
            
            powerMaps{i} = temp; 

        end
        
        fileList = dir(fullfile(savedir, sprintf('*%s*locationMap.mat', direction)));

        for ii=1:length(fileList)
            fname=strcat(savedir,fileList(ii).name);
            temp = cell2mat(struct2cell(load(fname)));
            
            %median filter
            filterRadius = 3; %note: must be odd
            nanMask = isnan(temp);
            temp_padded = padarray(temp, floor([filterRadius/2 filterRadius/2]), 'both'); % Pad image
            temp_paddedCol = im2col(temp_padded, [filterRadius filterRadius], 'sliding'); % Transform into columns
            temp_median = reshape(median(temp_paddedCol, 1, 'omitnan'), size(temp,1), size(temp,2)); % Find median of each column and reshape back
            temp_median(nanMask) = nan; % Set nan elements back
            
            locationMaps{i} = temp_median; 
        end
end

saveVariableList = {mean(cat(3, powerMaps{3:4}), 3); %azimuth_powerMap
                    mean(cat(3, locationMaps{3:4}), 3); %azimuth_locationMap
                    mean(cat(3, powerMaps{1:2}), 3); %altitude_powerMap
                    mean(cat(3, locationMaps{1:2}), 3)}; %altitude_locationMap

savefnameList = [strcat(savedir, "azimuth_powerMap.tif"), 
                 strcat(savedir, "azimuth_locationMap.tif"), 
                 strcat(savedir, "altitude_powerMap.tif"), 
                 strcat(savedir, "altitude_locationMap.tif")];

tagstruct.ImageLength = size(azimuth_powerMap,1);
tagstruct.ImageWidth = size(azimuth_powerMap,2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 64;
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP; 
tagstruct.SamplesPerPixel = 1;
tagstruct.Compression = Tiff.Compression.None;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

for i = 1:length(savefnameList)
        saveVar = cell2mat(saveVariableList(i));
        savefName = savefnameList(i);
        
        t = Tiff(char(savefName),'w');
        setTag(t,tagstruct);
        
        write(t, saveVar);
        close(t);
end



figure; imagesc(azimuth_locationMap, 'AlphaData', ~isnan(azimuth_locationMap)); 
set(gca,'color',[1 1 1]);
colormap(jet); colorbar; axis square
title('Azimuth location map','Interpreter','none');

figure; imagesc(azimuth_powerMap, 'AlphaData', ~isnan(azimuth_powerMap)); 
set(gca,'color',[1 1 1]);
colormap(jet); colorbar; axis square
title('Azimuth power map','Interpreter','none');

figure; imagesc(altitude_locationMap, 'AlphaData', ~isnan(altitude_locationMap)); 
set(gca,'color',[1 1 1]);
colormap(jet); colorbar; axis square
title('Altitude location map','Interpreter','none');

figure; imagesc(altitude_powerMap, 'AlphaData', ~isnan(altitude_powerMap)); 
set(gca,'color',[1 1 1]);
colormap(jet); colorbar; axis square
title('Altitude power map','Interpreter','none');
