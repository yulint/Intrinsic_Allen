clear all
%--------------------------------------------------------------------------
%% Define directories and files; change this section for every experiment

rootdir  = 'C:\Users\2018_Group_a\Desktop\yulin\intrinsic\';
fname = '20190207_172114';  
direction = 'L2R' ; %{'B2U', 'U2B', 'L2R' 'R2L'}

fname_BV = '';
fname_CAM = '';
fname_trigger = '';

savedir =  fullfile(rootdir, 'analysis_Allen\'); 
cursavedir = [savedir,fname,'\'];
if ~isdir(cursavedir), mkdir(cursavedir); end
%%
if isempty(fname_CAM)
    fname_CAM = [rootdir,fname,'trgfrm.bin'];
    clear CAM;

    CAM.fname = fname_CAM;
    CAM.fsize=dir(fname_CAM);
    fid = fopen(fname_CAM,'r');
    % framenumber uint32
    % baseline period over uint8
    % on period over uint8
    % timestamp double
    % interframe interval double

    CAM.frmsiz = (32/8+32/8+8/8+8/8+64/8+64/8+32/8); % size of single frame
    CAM.nfrm = (CAM.fsize.bytes)./CAM.frmsiz;

    % read in data: read all values of every parameter in at once
    fseek(fid,0,'bof'); % go to beginning of file
    CAM.imaqfcount       = fread(fid,CAM.nfrm,'*uint32',CAM.frmsiz-32/8,'ieee-be');

    fseek(fid,32/8,'bof');
    CAM.buf       = fread(fid,CAM.nfrm,'*int32',CAM.frmsiz-32/8,'ieee-be');

    fseek(fid,32/8+32/8,'bof'); % [double(CAM.i),double(CAM.buf+1),double(CAM.i)-double(CAM.buf+1)] isequal(double(CAM.i),double(CAM.buf+1)), find(diff(CAM.i~=0))
    CAM.imagingBslOver = fread(fid,CAM.nfrm,'*uint8',CAM.frmsiz-8/8,'ieee-be');

    fseek(fid,32/8+32/8+8/8,'bof');
    CAM.imagingStartEnd   = fread(fid,CAM.nfrm,'*uint8',CAM.frmsiz-8/8,'ieee-be');

    fseek(fid,32/8+32/8+8/8+8/8,'bof');
    CAM.timestamp= fread(fid,CAM.nfrm,'*double',CAM.frmsiz-64/8,'ieee-be');

    fseek(fid,32/8+32/8+8/8+8/8+64/8,'bof');
    CAM.timeifi  = fread(fid,CAM.nfrm,'*double',CAM.frmsiz-64/8,'ieee-be');

    fseek(fid,32/8+32/8+8/8+8/8+64/8+64/8,'bof');
    CAM.i  = fread(fid,CAM.nfrm,'*uint32',CAM.frmsiz-32/8,'ieee-be');        

    fclose(fid);

    save(sprintf('%s%s_CAM',cursavedir,fname),'CAM');

else
    CAM = load(sprintf('%s%s_CAM',cursavedir,fname'));
end 

% figure;plot(CAM.timestamp);
% figure;plot(CAM.timeifi);

%extract trigger info
if isempty(fname_trigger)
    fname_trigger = [rootdir,fname,'trgstm.bin'];
    clear TTL;

    TTL.fname = fname_trigger;
    TTL.fsize=dir(fname_trigger);
    fid = fopen(fname_trigger,'r');
    
    % trigid: 0,1 = sync pulses, 2 = trigger pulse
    % trgval: 1 for trigger pulse
    % i: frame when lines are high

    TTL.frmsiz = (32/8+8/8+8/8); % size of single frame
    TTL.nfrm = (TTL.fsize.bytes)./TTL.frmsiz;

    fseek(fid,0,'bof'); 
    TTL.i       = fread(fid,TTL.nfrm,'*uint32',TTL.frmsiz-32/8,'ieee-be');
    
    fseek(fid,32/8,'bof');
    TTL.trgid    = fread(fid,TTL.nfrm,'*uint8',TTL.frmsiz-8/8,'ieee-be');
    
    fseek(fid,32/8+8/8,'bof'); % [double(STM.i),double(STM.buf+1),double(STM.i)-double(STM.buf+1)] isequal(double(STM.i),double(STM.buf+1))
    TTL.trgval = fread(fid,TTL.nfrm,'*uint8',TTL.frmsiz-8/8,'ieee-be');
    
    fclose(fid);
    save(sprintf('%s%s_TTL',cursavedir,fname),'TTL');
else
    TTL = load(sprintf('%s%s_TTL',cursavedir,fname));
end 

% check TTL and imaging stamps

[pks,locs] = findpeaks(double(CAM.imagingStartEnd),'minpeakheight',0.1);

figure;
hold on;plot(CAM.i,CAM.imagingBslOver,'g')
hold on;plot(CAM.i,CAM.imagingStartEnd,'b');
hold on;plot(CAM.i(locs),0.5,'rx');
set(gca,'YLim',[-0.2 1.2]);
hold on;plot(TTL.i(TTL.trgid==0),0.25,'go')
hold on;plot(TTL.i(TTL.trgid==1),0.25,'gx')
hold on;plot(TTL.i(TTL.trgid==2),double(TTL.trgval(TTL.trgid==2))./10,'gs')
%%

%temporary version: TTLs recorded at end of baseline and end of cycle

fprintf('%s:%s: reading image\n',fname,datestr(now));

% read in image ------------------------------------------------------
fname_img=[rootdir,fname,'img.bin'];

clear IMG;
IMG.fname = fname_img;
IMG.fsize=dir(fname_img);

img_fid = fopen(fname_img,'r');

IMG.dim=fread(img_fid,2,'int32',0,'ieee-be');
IMG.avgnfrm=fread(img_fid,1,'int32',0,'ieee-be');

IMG.hdrsiz = 3*(32/8); % size of header
IMG.frminf = 32/8+32/8; % size of info written for each frame
IMG.frmsiz = IMG.frminf+( prod(IMG.dim)*16/8 ); % size of single frame
IMG.nfrm = (IMG.fsize.bytes-IMG.hdrsiz) ./IMG.frmsiz;

fseek(img_fid,IMG.hdrsiz,'bof');
IMG.i  = fread(img_fid,IMG.nfrm,'*uint32',IMG.frmsiz-32/8,'ieee-be');
fseek(img_fid,IMG.hdrsiz+32/8,'bof');
IMG.nf = fread(img_fid,IMG.nfrm,'*int32',IMG.frmsiz-32/8,'ieee-be');

% find peristimulus times  -------------------------------------------
clear TRIAL;
TRIAL.stimulusStarts = TTL.i(TTL.trgid==0);
TRIAL.cycleEnds = TTL.i(TTL.trgid==1);

%units here: number of CAM/TTL i-s  
TRIAL.avgTrialLength = double(median(TRIAL.cycleEnds(2:2:end) - TRIAL.cycleEnds(1:2:end))); % avg # i-s in each trial
TRIAL.avgPreGap = double(median(TRIAL.stimulusStarts(2:end)-TRIAL.cycleEnds(1:end-1)));
TRIAL.avgPostGap = 3 * 1./median(diff(CAM.timestamp));

%units here: number of IMG frames
%find # baseline frames 
% TRIAL.preBaselineTime = 3;
% PSTH.frmrate = 24;
 IMG.avgfrmTime = median(diff(CAM.timestamp(IMG.i))); % time elapsed in average i step of IMG
 PSTH.frmrate = 1./IMG.avgfrmTime; %error occurs when frames droped? (20.10.2014.)
% PSTH.frmpre  = round(TRIAL.preBaselineTime*PSTH.frmrate);

% extract PSTH  -------------------------------------------

fprintf('%s:%s: extracting PSTH\n',fname,datestr(now));

% create cell array of movies ([X Y frames] arrays), one per cycle
PAllTrials = cell(1,length(TRIAL.stimulusStarts));
PAllTrials_dFoverF = cell(1,length(TRIAL.stimulusStarts));

for trialNo=1:length(TRIAL.stimulusStarts) 
    fprintf('trial number: %d\n',trialNo);
    currentTrialStart = TRIAL.stimulusStarts(trialNo) - TRIAL.avgPreGap;
    currentTrial = [currentTrialStart,currentTrialStart+TRIAL.avgTrialLength]; 
    currentBaselineEnd = TRIAL.stimulusStarts(trialNo);

    % find frame with starttime closest to trigger
    [~,frIx1]=min(abs(double(IMG.i)-double(currentTrial(1))));
    [~,frIx2]=min(abs(double(IMG.i)-double(currentTrial(2))));
    frmrng = frIx1:frIx2;

    %find baseline frames, relative to start of trial
    [~,frIxBsl]=min(abs(double(IMG.i)-double(currentBaselineEnd)));
    bslRng = [1:frIxBsl-frIx1]; 

    if ~any(frmrng<1|frmrng>IMG.nfrm)
        nfrm = length(frmrng);

        skip=IMG.hdrsiz+(frmrng(1)-1).*IMG.frmsiz + IMG.frminf; % header size + preceding frames + frm info of first frame
        fseek(img_fid,skip,'bof'); % go to beginning of first frame
        I=fread(img_fid,nfrm*prod(IMG.dim),sprintf('%d*int16',prod(IMG.dim)),IMG.frminf,'ieee-be');

        TI = reshape(I,[IMG.dim(1),IMG.dim(2),nfrm]); % size(TI)

        PAllTrials{trialNo}  = TI;

        baselineF =mean(TI(:,:,bslRng),3); 
        baselineF =repmat(baselineF,[1,1,size(TI,3)]);
        dF=(TI-baselineF)./baselineF;
        PAllTrials_dFoverF{trialNo}=dF;

    else
        error('Problem Matching the Frames')
    end

end

fclose(img_fid);

%% average movie across trials, matching # frames to the trial with fewest
% frames

min_nfr = min(cellfun(@(x) size(x,3) , PAllTrials_dFoverF)); % min_nfr: minimum num frames across the movies

%PAveraged=NaN([length(PAllTrials),size(PAllTrials{trialNo},1),size(PAllTrials{trialNo},2),min_nfr]);
PAveraged_dFoverF=NaN([length(PAllTrials_dFoverF),size(PAllTrials_dFoverF{trialNo},1),size(PAllTrials_dFoverF{trialNo},2),min_nfr]);
for trialNo=1:length(PAllTrials_dFoverF)
    %PAveraged(trialNo,:,:,:)=PAllTrials{trialNo}(:,:,(1:min_nfr));
    PAveraged_dFoverF(trialNo,:,:,:)=PAllTrials_dFoverF{trialNo}(:,:,(1:min_nfr));
end
%PAveraged=squeeze(nanmean(PAveraged,1));
PAveraged_dFoverF=squeeze(mean(PAveraged_dFoverF,1), 'omitnan');


%clear PAllTrials % large cell array of movies, not used anymore
%clear PAllTrials_dFoverF
%% Fourier analysis
fprintf('%s:%s: fft\n',fname,datestr(now));

% fft across time; dimensions of frame_fft = rows of pixels, columns of
% pixels, fft frequency
frame_fft = fft(PAveraged_dFoverF,[],3);
freq_range = (0:min_nfr-1)*(PSTH.frmrate/min_nfr);

%find first harmonic frequency, which should also be the frequency of the
%visual stimulus sweep, and which should also be the first non-zero
%frequency in fft (since movie has been cut to length of each sweep) 
visualSweepFreq_img = 1/ (min_nfr/PSTH.frmrate); %based on number of imaged frames
visualSweepFreq_ttl = 1/ double(TRIAL.avgTrialLength * median(diff(CAM.timestamp))); % based on number of 'i's between TTL sync pulses
firstNonZeroFreq_fft = freq_range(2);

fprintf('first non-zero freq in fft = %d \n check that similar to visual sweep freq based on # imaged frames: %d, \nbased on sync pulses: %d \n ',visualSweepFreq_img,visualSweepFreq_ttl,firstNonZeroFreq_fft);

fft_firstHarmonicFreq = frame_fft(:,:,2);

% power map
powerMovie = (abs(frame_fft)*2)/ min_nfr; %normalise by number of frames
powerMap = abs(powerMovie(:,:,2)');
%powerMap = abs(fft_firstHarmonicFreq'); %transpose to match back to BV image

%phase map
% phaseMovie = angle(frame_fft);
% phaseMap = -1 *  phaseMovie(:,:,2)';
% phaseMap = mod(phaseMap, 2 *pi);
phaseMap = 180/pi* mod(angle(fft_firstHarmonicFreq'), 2*pi);

%sanity check
%figure; h = heatmap(phaseMap, "GridVisible", "Off", "Colormap", jet);
fig_pwr = figure; imagesc(powerMap); colormap(gray); colorbar; axis square
title(sprintf('Power map_%s_%s', fname, direction),'Interpreter','none');
powerMap_fname = sprintf('%s%s_%s_powerMap',cursavedir,fname,direction);
print(fig_pwr, powerMap_fname, '-dtiff');

%%
%convert phase map to location map
%relative location in screen, 0 = bottom/ left of screen

TRIAL.visualStimStartPhase = (TRIAL.avgPreGap)/TRIAL.avgTrialLength * 360;
TRIAL.visualStimEndPhase = (TRIAL.avgTrialLength-TRIAL.avgPostGap)/TRIAL.avgTrialLength * 360;

locationMap = phaseMap;
locationMap(locationMap < TRIAL.visualStimStartPhase | locationMap > TRIAL.visualStimEndPhase) = NaN;

for index = find(locationMap > TRIAL.visualStimStartPhase & locationMap < TRIAL.visualStimEndPhase)
    if strcmp(direction, 'B2U') || strcmp(direction, 'L2R')
        locationMap(index) = (locationMap(index)- TRIAL.visualStimStartPhase)/ (TRIAL.avgTrialLength-TRIAL.avgPreGap-TRIAL.avgPostGap);
    elseif strcmp(direction, 'U2B') || strcmp(direction, 'R2L')
        locationMap(index) = 1 - (locationMap(index)- TRIAL.visualStimStartPhase)/ (TRIAL.avgTrialLength-TRIAL.avgPreGap-TRIAL.avgPostGap);
    end    
end

% figure; h_location = heatmap(locationMap, "GridVisible", "Off", "Colormap", jet);
% h_location.Title = sprintf('Location map_%s_%s', fname, direction); 

fig_locMap = figure; imagesc(locationMap, 'AlphaData', ~isnan(locationMap)); 
set(gca,'color',[1 1 1]);
colormap(jet); colorbar; axis square
title(sprintf('Location map_%s_%s', fname, direction),'Interpreter','none');
locationMap_fname = sprintf('%s%s_%s_locationMap',cursavedir,fname,direction);
print(fig_locMap, locationMap_fname, '-dtiff');

%%
save(sprintf('%s%s_%s_fft',savedir,fname,direction),'fft_firstHarmonicFreq')
save(sprintf('%s%s_%s_powerMap',savedir,fname,direction),'powerMap')
save(sprintf('%s%s_%s_locationMap',savedir,fname,direction),'locationMap')

