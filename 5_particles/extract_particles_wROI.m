% Read Data
clear, clc, close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IM_number = 18;

fitting_dist        = '/Users/christian/Documents/Arbeit/MatLab/SMLM_tutorial/localizations';

Locpath1            = ['/Volumes/Christian-Sieben/data_HTP/2019-04-02_Inflammasome_3D_ASC_speck_Ab/locResults_Feng/THP_1_' num2str(IM_number) '_1'];
locName1            = ['THP_1_' num2str(IM_number) '_1_MMStack_1_Localizations_Z'];

WFpath              = ['/Volumes/Christian-Sieben/data_HTP/2019-04-02_Inflammasome_3D_ASC_speck_Ab/THP_1_WF' num2str(IM_number)];
WFname              = ['THP_1_WF' num2str(IM_number) '_MMStack_Pos0.ome.tif'];

resultsFolder       = ['/Volumes/Transcend/Inflammasome/2020-03-02_CS_Inflammasome_3D_Ab'];
resultsFile         = '2019-04-02_CS_Inflammasome_3D_Ab_Fang_Particles_w_clusters.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pixelsize 	= 106; % nm

fprintf('\n -- Path and File information loaded --\n')
cd(Locpath1);

% Load data

cd(Locpath1);
locs=dlmread([locName1 '.csv'],',',1,0);

fprintf('\n -- Data loaded --\n')
% Load header

cd(Locpath1);
file    = fopen([locName1 '.csv']);
line    = fgetl(file);
h       = regexp( line, ',', 'split' );

% xCol                = strmatch('"x [nm]"',h);
% yCol                = strmatch('"y [nm]"',h);
% zCol                = strmatch('"z [nm]"',h);
% photonsCol          = strmatch('"intensity [photons]"',h);
% framesCol           = strmatch('"frame"',h);
% LLCol               = strmatch('"loglikelihood"',h);

xCol                = strmatch('x [nm]',h);
yCol                = strmatch('y [nm]',h);
zCol                = strmatch('z [nm]',h);
photonsCol          = strmatch('intensity [photon]',h);
framesCol           = strmatch('frame',h);
LLCol               = strmatch('loglikelihood',h);
BG_col              = strmatch('background',h);

% xCol                = strmatch('x_nm',h);
% yCol                = strmatch('y_nm',h);
% xCol_px             = strmatch('x_pix',h);
% yCol_px             = strmatch('y_pix',h);
% zCol                = strmatch('z_nm',h);
% photonsCol          = strmatch('photons',h);
% framesCol           = strmatch('frame',h);
% LLCol               = strmatch('logLikelyhood',h);
% BG_col              = strmatch('background',h);

% Load WF image
cd(WFpath)
WF_image    = imread(WFname);

% Show histograms of filter parameters

fprintf('\n -- Plotting histograms --\n')
close all

figure('Position',[50 300 500 500])
subplot(2,2,1);
bins = -1000:50:0;
h1 = hist(locs(:,LLCol),bins);
bar(bins,h1/sum(h1));
xlabel('LL')
title(['Median = ' num2str(median(locs(:,LLCol)))]);
axis([-1000 0 0 max(h1/sum(h1))*1.2])


subplot(2,2,2);
bins = 0:200:10000;
h1 = hist(locs(:,photonsCol),bins);
bar(bins,h1/sum(h1));
xlabel('Photons');
title(['Median = ' num2str(median(locs(:,photonsCol)))]);
axis([0 1e4 0 max(h1/sum(h1))*1.2])

if isempty(zCol)==1;
    
subplot(2,2,3);
title('No Z column found');
    
else

subplot(2,2,3);
bins = -2000:200:2000;
h1 = hist(locs(:,zCol),bins);
bar(bins,h1/sum(h1));
xlabel('Z position');
title(['Median = ' num2str(median(locs(:,zCol)))]);
axis([-2000  2000 0 max(h1/sum(h1))*1.2])

end

if isempty(BG_col)==1;
    
subplot(2,2,4);
title('No BG column found');
    
else

subplot(2,2,4);
bins = 0:5:150;
h1 = hist(locs(:,BG_col),bins);
bar(bins,h1/sum(h1));
xlabel('BG');
title(['Median = ' num2str(median(locs(:,BG_col)))]);
axis([0 150 0 max(h1/sum(h1))*1.2])

end

figure('Position',[550 300 500 500])
imagesc(flipud(WF_image),[0 10*median(WF_image(:))]);hold on;
axis square

%% Drift Correction

% Generate Coords variable as input for RCC
% Remove NaN and Inf

clc
cd([fitting_dist '/RCC_drift_correct']);

% coords(:,1) = locs(:,xCol_px);
% coords(:,2) = locs(:,yCol_px);
% coords(:,3) = locs(:,framesCol);

coords(:,1) = locs(:,xCol)/pixelsize;
coords(:,2) = locs(:,yCol)/pixelsize;
coords(:,3) = locs(:,framesCol);

% Remove NaN and Inf
temp = coords;
clear coords
coords = temp( ~any( isnan( temp(:,1) ) | isinf( temp(:,1) ), 2 ),: );

fprintf('\n -- Ready for DC -- \n')

% Select region to correct

pxlsize = 1;

heigth  = round((max(coords(:,2))-min(coords(:,2)))/pxlsize);
width   = round((max(coords(:,1))-min(coords(:,1)))/pxlsize);
im      = hist3([coords(:,1),coords(:,2)],[width heigth]); % heigth x width

% Select rectangles

rect = []; 

figure('Position',[100 200 600 600])
f = imagesc(imrotate(im,90),[0 20]);
colormap('parula'); colorbar;

rect = getrect;

fprintf('\n -- ROI selected --\n')

xmin = min(coords(:,1))+ rect(1,1)*pxlsize;
ymin = max(coords(:,2))- rect(1,2)*pxlsize - (rect(1,4)*pxlsize) ;
xmax = xmin + 512;
ymax = ymin + 512;

target      = find(coords(:,1)>xmin & coords(:,1)<xmax & coords(:,2)>ymin & coords(:,2)<ymax);
coords_ROI  = coords(target,1:end);

% Show cropped region

heigth  = round((max(coords_ROI(:,2))-min(coords_ROI(:,2)))/pxlsize);
width   = round((max(coords_ROI(:,1))-min(coords_ROI(:,1)))/pxlsize);
im      = hist3([coords_ROI(:,1),coords_ROI(:,2)],[width heigth]); % heigth x width

figure('Position',[100 200 600 600])
imagesc(imrotate(im,90),[0 20]);
colormap('parula'); colorbar;

coords_ROI(:,1) = coords_ROI(:,1)-min(coords_ROI(:,1));
coords_ROI(:,2) = coords_ROI(:,2)-min(coords_ROI(:,2));


% Drift correct
tic
close all
% Input:    coords:             localization coordinates [x y frame], 
%           segpara:            segmentation parameters (time wimdow, frame)
%           imsize:             image size (pixel)
%           pixelsize:          camera pixel size (nm)
%           binsize:            spatial bin pixel size (nm)
%           rmax:               error threshold for re-calculate the drift (pixel)
% Output:   coordscorr:         localization coordinates after correction [xc yc] 
%           finaldrift:         drift curve (save A and b matrix for other analysis)

%finaldrift = RCC_TS(filepath, 1000, 256, 160, 30, 0.2);

segpara     = max(locs(:,framesCol))/5; 
imsize      = 512;
pixelsize 	= 106;
binsize     = 30;

[coordscorr, finaldrift] = RCC(coords_ROI, segpara, imsize, pixelsize, binsize, 0.2);
% [coordscorr, finaldrift] = DCC(coords, segpara, imsize, pixelsize, binsize);
                          
clc

display(['DC finished in ' round(num2str(toc/60)) ' min']);

% Plot Drift curves

figure('Position',[100 100 900 400])
subplot(2,1,1)
plot(finaldrift(:,1))
title('x Drift')
subplot(2,1,2)
plot(finaldrift(:,2))
title('y Drift')

%% Apply correction to coords

deltaXY      = [];
deltaXY(:,1) = coords(:,3); % frames
deltaXY(:,2) = finaldrift(deltaXY(:,1),1); % x drift
deltaXY(:,3) = finaldrift(deltaXY(:,1),2); % y drift

coordsDC      = [];
coordsDC(:,1) = coords(:,1)-deltaXY(:,2);
coordsDC(:,2) = coords(:,2)-deltaXY(:,3);

display('Coord corrected');
   
% Generate DC variable

locs_DC = locs;
% locs_DC(:,xCol_px)  = coordsDC(:,1); % x pxl
% locs_DC(:,yCol_px)  = coordsDC(:,2); % y pxl
% locs_DC(:,xCol) = coordsDC(:,1)*pixelsize; % x nm
% locs_DC(:,yCol) = coordsDC(:,2)*pixelsize; % y nm
% 
locs_DC(:,xCol) = coordsDC(:,1)*pixelsize; % x nm
locs_DC(:,yCol) = coordsDC(:,2)*pixelsize; % y nm


close all


%% Filter localizations

minFrame            = 1000;
maxFrame            = max(locs(:,framesCol));
MinPhotons          = 1000;
MaxPhotons          = 100000;
MaxLL               = -500;
z_min               = -800;
z_max               = 500;
BG_max              = 100; 
        
filter              = [];
filter              = find(locs_DC(:,photonsCol) > MinPhotons & locs_DC(:,photonsCol) < MaxPhotons & locs_DC(:,LLCol) > MaxLL ... 
                         & locs_DC(:,framesCol) > minFrame & locs_DC(:,4) > z_min & locs_DC(:,4) < z_max & locs_DC(:,BG_col) < BG_max ...
                         & locs_DC(:,framesCol) < maxFrame);
                     
locsFilt            = locs_DC(filter,1:end);

clc
display(['Localizations filtered (' num2str(length(locsFilt)/length(locs_DC)) ' left)']);


%% Out 1a. Extract particle locs in ROI

cd('/Users/christian/Documents/Arbeit/MatLab/SMLM_tutorial/Particles')

figure('Position',[50 500 400 400])
imagesc(flipud(WF_image),[0 10*median(WF_image(:))]);hold on;

[locs_ROI] = manualROI(locs_DC,xCol,yCol,100);

ParticlesTemp = {};
ParticlesTemp{1,1} = locs_ROI;

fprintf('\n -- Out 1 Particle locs extracted -- \n')

%% Out 1b. Save ParticlesTemp

% First time
% ParticlesTemp = {};
% ParticlesTemp = cat(1, ParticlesTemp, locs_ROI);

% After first time 
% cd(resultsFolder);load(resultsFile);
% ParticlesTemp = cat(1, ParticlesTemp, locs_ROI);

%% Out 2. Crop WF image

% fig = imshow(WF_image)
figure
imagesc(WF_image,[0 10*median(WF_image(:))]);hold on;
[xi,yi] = getpts

croppedWF = imcrop(WF_image,[xi-25 yi-25 50 50]);
% imshow(croppedWF)
figure
imagesc(croppedWF,[0 10*median(croppedWF(:))]);hold on;

ParticlesTemp{size(ParticlesTemp, 1),2} = croppedWF;

fprintf('\n -- Out 2 WF image cropped -- \n')

%% Out 3. DBSCAN with ROI

cd('/Users/christian/Documents/Arbeit/MatLab/SMLM_tutorial/localizations');

[subset, LpA, CpA, LpC] = DBSCAN_with_ROI(locs_DC,xCol, yCol, 200, 5, 50);
 
% Clusters per Area
% Locs per area

ParticlesTemp{size(ParticlesTemp, 1),3} = subset;
ParticlesTemp{size(ParticlesTemp, 1),4} = LpC; % locs per cluster
ParticlesTemp{size(ParticlesTemp, 1),5} = CpA; % clusters per area/ µm2
ParticlesTemp{size(ParticlesTemp, 1),6} = IM_number; % Image ID

fprintf('\n -- 3 Cytoplasm clustered -- \n')

%% Save results
close all
cd(resultsFolder);load(resultsFile);
% Particles = ParticlesTemp;
Particles = cat(1, Particles, ParticlesTemp);
save(resultsFile,'Particles');

clc
display('Particle file updated');
