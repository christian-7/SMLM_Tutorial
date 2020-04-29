%% Open localization file

clear, clc, close all

fitting_dist = '/Users/christian/Documents/Arbeit/MatLab/SMLM_tutorial/localizations';                                  % define path to fitting distribution

path        = '/Volumes/Transcend/data_HTP/2019-05-09_CS_3D_Microtubules';
savepath    = '/Volumes/Transcend/data_HTP/2019-05-09_CS_3D_Microtubules/analysis';

cd(path);

filename = 'Cos7_MT_A647_4_1_locs';
locs     = dlmread([filename '.csv'],',',1,0);

file          = fopen([filename '.csv']);
line          = fgetl(file);
header        = regexp( line, ',', 'split' );

pixelsize 	= 110; % nm

xCol            = strmatch('x_nm',header);
yCol            = strmatch('y_nm',header);
z_nm            = strmatch('z_nm',header);
framesCol       = strmatch('frame',header);
LLCol           = strmatch('logLikelyhood',header);
photonsCol      = strmatch('photons',header);
xCol_px         = strmatch('x_pix',header);
yCol_px         = strmatch('y_pix',header);
BG_col          = strmatch('crlb_background',header);

fprintf('\n -- Loading file ... \n')

fprintf('\n -- Data Loaded -- \n')

%% Show histograms of filter parameters

close all

figure('Position',[400 300 700 700])
subplot(2,2,1);
bins = -1000:50:0;
h1 = hist(locs(:,LLCol),bins);
bar(bins,h1/sum(h1));
xlabel('LL')
title(['Median = ' num2str(median(locs(:,LLCol)))]);
axis([-1000 0 0 max(h1/sum(h1))])


subplot(2,2,2);
bins = 0:200:10000;
h1 = hist(locs(:,photonsCol),bins);
bar(bins,h1/sum(h1));
xlabel('Photons');
title(['Median = ' num2str(median(locs(:,photonsCol)))]);
axis([0 1e4 0 max(h1/sum(h1))])

if isempty(z_nm)==1;
    
subplot(2,2,3);
title('No Z column found');
    
else

subplot(2,2,3);
bins = -2000:200:2000;
h1 = hist(locs(:,z_nm),bins);
bar(bins,h1/sum(h1));
xlabel('Z position');
title(['Median = ' num2str(median(locs(:,z_nm)))]);
axis([-2000  2000 0 max(h1/sum(h1))])

end

subplot(2,2,4);
bins = 0:0.2:5;
h1 = hist(locs(:,BG_col),bins);
bar(bins,h1/sum(h1));
xlabel('BG');
title(['Median = ' num2str(median(locs(:,BG_col)))]);
axis([-2 10 0 max(h1/sum(h1))])


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%START%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Drift Correction%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate Coords variable as input for RCC
% Remove NaN and Inf

clc
cd([fitting_dist '/RCC_drift_correct']);

coords(:,1) = locs(:,xCol_px);
coords(:,2) = locs(:,yCol_px);
coords(:,3) = locs(:,framesCol);

% Remove NaN and Inf
temp = coords;
clear coords
coords = temp( ~any( isnan( temp(:,1) ) | isinf( temp(:,1) ), 2 ),: );

fprintf('\n -- Ready for DC -- \n')


%% Select region to correct

pxlsize = 1;

heigth  = round((max(coords(:,2))-min(coords(:,2)))/pxlsize);
width   = round((max(coords(:,1))-min(coords(:,1)))/pxlsize);
im      = hist3([coords(:,1),coords(:,2)],[width heigth]); % heigth x width

% Select rectangles

rect = []; 

figure('Position',[100 200 600 600])
f = imagesc(imrotate(im,90),[0 2e3]);
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
imagesc(imrotate(im,90),[0 2e3]);
colormap('parula'); colorbar;

coords_ROI(:,1) = coords_ROI(:,1)-min(coords_ROI(:,1));
coords_ROI(:,2) = coords_ROI(:,2)-min(coords_ROI(:,2));


%% Drift correct
tic

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

%% Plot Drift curves

figure('Position',[100 100 900 400])
subplot(2,1,1)
plot(finaldrift(:,1))
title('x Drift')
subplot(2,1,2)
plot(finaldrift(:,2))
title('y Drift')

%% Apply correction to coords

deltaXY = [];
deltaXY(:,1) = coords(:,3); % frames
deltaXY(:,2) = finaldrift(deltaXY(:,1),1); % x drift
deltaXY(:,3) = finaldrift(deltaXY(:,1),2); % y drift

coordsDC      = [];
coordsDC(:,1) = coords(:,1)-deltaXY(:,2);
coordsDC(:,2) = coords(:,2)-deltaXY(:,3);

display('Coord corrected');
   
% Generate DC variable

locs_DC = locs;
locs_DC(:,xCol_px)  = coordsDC(:,1); % x pxl
locs_DC(:,xCol_px)  = coordsDC(:,2); % y pxl
locs_DC(:,xCol) = coordsDC(:,1)*pixelsize; % x nm
locs_DC(:,yCol) = coordsDC(:,2)*pixelsize; % y nm

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Drift Correction%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Filter localizations

minFrame            = 1000;
maxFrame            = max(locs(:,framesCol));
MinPhotons          = 2000;
MaxPhotons          = 100000;
MaxLL               = -600;
z_min               = -500;
z_max               = 500;
BG_max              = 7; 
        
filter              = [];
filter              = find(locs_DC(:,photonsCol) > MinPhotons & locs_DC(:,photonsCol) < MaxPhotons & locs_DC(:,LLCol) > MaxLL ... 
                         & locs_DC(:,framesCol) > minFrame & locs_DC(:,4) > z_min & locs_DC(:,4) < z_max & locs_DC(:,BG_col) < BG_max ...
                         & locs_DC(:,framesCol) < maxFrame);
                     
locsFilt            = locs_DC(filter,1:end);

clc
display(['Localizations filtered (' num2str(length(locsFilt)/length(locs_DC)) ' left)']);

%% Render localizations

% locsFilt = coordsDC; % xCol = 1; yCol = 2;

pxlsize = 20; % nm

heigth = round((max(locsFilt(:,yCol)) - min((locsFilt(:,yCol))))/pxlsize);
width  = round((max(locsFilt(:,xCol)) - min((locsFilt(:,xCol))))/pxlsize);
        
rendered = hist3([locsFilt(:,yCol),locsFilt(:,xCol)],[heigth width]);

figure 
%imagesc(imadjust(imgaussfilt(rendered,1)));
imshow(imgaussfilt(rendered,1),[0.5 1]);
colormap hot


% figure 
% scatter3(locsFilt(:,xCol),locsFilt(:,yCol),locsFilt(:,z_nm),1);

%% Generate output for ThunderSTORM

tic;
cd(savepath);

locsTS = [];
locsTS(:,1) = locsFilt(:,1);                % frames
locsTS(:,2) = locsFilt(:,13);               % x nm
locsTS(:,3) = locsFilt(:,14);               % y nm
locsTS(:,4) = locsFilt(:,4);                % z nm
locsTS(:,5) = locsFilt(:,6);                % photons
locsTS(:,6) = locsFilt(:,12);               % LL

NameCorrected = [filename '_DC_forTS.csv'];

fileID = fopen(NameCorrected,'w');
fprintf(fileID,[['"frame","x [nm]","y [nm]","z [nm]","intensity [photons]","loglikelihood"'] ' \n']);
dlmwrite(NameCorrected,locsTS,'-append');
fclose('all');

clc;

clc
display(['File saved in ' num2str(toc) ' min']);


%% Save image

cd(savepath);

I32 = [];
I32 = uint32(rendered);

name = [filename '_' num2str(pxlsize) '_nm_pxl.tiff'];

t = Tiff(name,'w');
tagstruct.ImageLength     = size(I32,1);
tagstruct.ImageWidth      = size(I32,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 32;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct)

t.write(I32);
t.close()

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Select ROI to export %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pxlsize = 50; locs_Filt_ROI = [];

heigth  = round((max(locsFilt(:,13))-min(locsFilt(:,13)))/pxlsize);
width   = round((max(locsFilt(:,14))-min(locsFilt(:,14)))/pxlsize);
im      = hist3([locsFilt(:,13),locsFilt(:,14)],[width heigth]); % heigth x width

% Select rectangles

rect = []; 

figure('Position',[100 200 600 600])
f = imagesc(imrotate(im,90), [0 1e2]); 
colormap('hot'); colorbar;

rect = getrect; % [xmin ymin width height].

fprintf('\n -- ROI selected --\n')

xmin = min(locsFilt(:,13))+ rect(1,1)*pxlsize;
ymin = max(locsFilt(:,14)) - rect(1,2)*pxlsize - (rect(1,4)*pxlsize) ;
xmax = xmin + (rect(1,3)* pxlsize);
ymax = ymin + (rect(1,4) * pxlsize);


target      = find(locsFilt(:,13)>xmin & locsFilt(:,13)<xmax & locsFilt(:,14)>ymin & locsFilt(:,14)<ymax);
locs_Filt_ROI  = locsFilt(target,1:end);

% Show cropped region

heigth  = round((max(locs_Filt_ROI(:,14))-min(locs_Filt_ROI(:,14)))/10);
width   = round((max(locs_Filt_ROI(:,13))-min(locs_Filt_ROI(:,13)))/10);
im      = hist3([locs_Filt_ROI(:,13),locs_Filt_ROI(:,14)],[width heigth]); % heigth x width

figure('Position',[100 200 600 600])
imagesc(imrotate(im,90),[0 5]);
colormap('hot'); colorbar


%% VISP format x,y,z,xPrec, yPrec, zPrec, Photons, Frame

cd(savepath);
locs_subset = [];

locs_subset (:,1) = locs_Filt_ROI(:,xCol);
locs_subset (:,2) = locs_Filt_ROI(:,yCol);
locs_subset (:,3) = locs_Filt_ROI(:,z_nm);
locs_subset (:,4) = 10;
locs_subset (:,5) = 10;
locs_subset (:,6) = 30;
locs_subset (:,7) = locs_Filt_ROI(:,photonsCol);
locs_subset (:,8) = locs_Filt_ROI(:,framesCol);

name_for_VISP = [filename, '_ROI_2_forVISP.txt'];

dlmwrite(name_for_VISP, locs_subset,'delimiter', '\t');

fprintf('\n -- Data saved in VISP format --\n')

%% Write for cake23
% https://www.cake23.de/pointcloud-loader/

cd(savepath);
locs_Flexi = [];

locs_Flexi(:,1) = locs_Filt_ROI(:,framesCol);
locs_Flexi(:,2) = locs_Filt_ROI(:,xCol);
locs_Flexi(:,3) = locs_Filt_ROI(:,yCol);
locs_Flexi(:,4) = locs_Filt_ROI(:,z_nm);
locs_Flexi(:,5) = 5;
locs_Flexi(:,6) = 10;

NameCorrected = [filename '_ROI_1_DC_cake23.txt'];

fileID = fopen(NameCorrected,'w');
fprintf(fileID,[['frame, x [nm], y [nm], z [nm],uncertainty_xy [nm], uncertainty_z [nm] '] ' \n']);
dlmwrite(NameCorrected,locs_Flexi,'-append','delimiter', '\t');
fclose('all');

% dlmwrite(NameCorrected,locs_Flexi,'delimiter', '\t');

%% Generate output for ThunderSTORM

tic;
cd(savepath);

locsTS = [];
locsTS(:,1) = locs_Filt_ROI(:,1);                % frames
locsTS(:,2) = locs_Filt_ROI(:,13);               % x nm
locsTS(:,3) = locs_Filt_ROI(:,14);               % y nm
locsTS(:,4) = locs_Filt_ROI(:,4);                % z nm
locsTS(:,5) = locs_Filt_ROI(:,5);                % photons
locsTS(:,6) = locs_Filt_ROI(:,12);               % LL

NameCorrected = [filename '_ROI_2_spline.csv'];

fileID = fopen(NameCorrected,'w');
fprintf(fileID,[['"frame","x [nm]","y [nm]","z [nm]","intensity [photons]","loglikelihood"'] ' \n']);
dlmwrite(NameCorrected,locsTS,'-append');
fclose('all');

clc;

clc
display(['File saved in ' num2str(toc) ' min']);

