%% sptPALM tracking and post-processing

% Author: Christian Sieben, EPFL 
% sieben.christian@gmail.com
% January 2020

% Content

% 1. Load loc file 
% 2. Filter and vizualize
% 3. Select ROI
% 4. Track using C&G  
% 5. Filter tracks
% 6. Export tracks

%% 1. Load Localizations

clear, close all, clc, clear

pathd = ('/Users/csi20/Documents/Transcend/data_PALM/2020-03-05_M1_mEos/');            % Data path
pathR = ('/Users/csi20/Documents/Transcend/data_PALM/2020-03-05_M1_mEos/analysis/');   % Results path

imageID = 4;

name_base = ['A549_M1mEosC1_447_gain300_30ms_00' num2str(imageID)];

pxl_size    = 160;  % nm
time_step   = 0.03; % in sec

%%%%%%%%%% Load File %%%%%%%%%%

cd(pathd);
locs        = dlmread([name_base '.csv'],',',1,0);
WF_image    = imread(['STD_' name_base '.tif']);                        % or STD image

file = fopen([name_base '.csv']);
line = fgetl(file);
h = regexp( line, ',', 'split' );

% Data from Matlab (Fangs, fitcspline)

xCol        = strmatch('x [nm]',h);
yCol        = strmatch('y [nm]',h);
frameCol    = strmatch('frame',h);
sigmaCol    = strmatch('sigma [nm]',h);
photonsCol  = strmatch('intensity [photon]',h);

% From Thunderstorm

xCol        = strmatch('"x [nm]"',h);
yCol        = strmatch('"y [nm]"',h);
frameCol    = strmatch('"frame"',h);
sigmaCol    = strmatch('"sigma [nm]"',h);
photonsCol  = strmatch('"intensity [photon]"',h);

fprintf('\n -- Data loaded --\n')

%% 2. Filter and generate PosList

close all

filtered    = []; density =[];
NNdist      = 400; % nm for the density plot
NoNind      = [];

%%%%% Filter


min_int     = 10;   % min photons
min_sigma   = 100;  % min sigma
max_sigma   = 200;  % max sigma


filter          = find(locs(:,photonsCol)>min_int & locs(:,sigmaCol)>min_sigma & locs(:,sigmaCol)<max_sigma);
filtered        = locs(filter,xCol);
filtered(:,2)   = locs(filter,yCol);
filtered(:,3)   = locs(filter,frameCol);

%%%%% Density

density = filtered(:,1:2);

idx = rangesearch(density,density,NNdist);
        
for i=1:length(idx);
NoNind(i,1)=length(idx{i,1});     % count the total number of neighbors for each point in dataset
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Position',[100 600 600 600],'Name',['STD image']);
imagesc(flipud(WF_image),[0 10*median(WF_image(:))]);hold on;
colormap(gray)
axis square;
xlabel('pxl');
ylabel('pxl');
set(gca,'FontSize',12)

figure('Position',[100 600 1200 400])
subplot(1,2,1)
scatter(filtered(:,1), filtered(:,2),1); hold on;
title([num2str(length(filtered)),' from ',num2str(length(locs)),' left after filtering -> ', num2str(length(filtered)/length(locs)*100),'  %',]);
axis([min(filtered(:,1)) max(filtered(:,1)) min(filtered(:,2)) max(filtered(:,2))])
xlabel('nm');
ylabel('nm');
set(gca,'FontSize',12)
box on; axis equal;

subplot(1,2,2)
scatter(density(:,1),density(:,2),5,NoNind,'o','filled');
colormap(jet)
title(['Localization density, NN Dist = ', num2str(NNdist) ' nm']);
axis([min(density(:,1)) max(density(:,1)) min(density(:,2)) max(density(:,2))])
xlabel('nm');
ylabel('nm');
whitebg(1,'k')
set(gca,'FontSize',12)
box on; axis equal;

% Create pos_list for track.m --> tracking in pxl/frames

pos_list=[];

pos_list(:,1)=filtered(:,1)/pxl_size;                   % in pxl
pos_list(:,2)=filtered(:,2)/pxl_size;                   % in pxl
pos_list(:,3)=filtered(:,3);                            % dt in frames

fprintf('\n -- Data filtered / PosList generated--\n')

%% 3. Select ROI
close all

pxlsize = 0.5;

heigth  = round((max(pos_list(:,2))-min(pos_list(:,2)))/pxlsize);
width   = round((max(pos_list(:,1))-min(pos_list(:,1)))/pxlsize);
im      = hist3([pos_list(:,1),pos_list(:,2)],[width heigth]); % heigth x width

% Select rectangles

rect = [];

figure('Position',[100 200 600 600],'Name','Select ROI for tracking')
imagesc(imrotate(im,90),[0 10]);
title('Please select ROI');
xlabel('x [pxl]');
ylabel('y [pxl]');
set(gca,'FontSize',12)
axis square

rect = getrect;

fprintf('\n -- ROI selected --\n')

% Plot fiducials and average curve rectangles
   
xmin = min(pos_list(:,1))+ rect(1)*pxlsize;
ymin = max(pos_list(:,2)) - rect(2)*pxlsize - (rect(4)*pxlsize) ;
xmax = xmin + (rect(3)* pxlsize);
ymax = ymin + rect(4) * pxlsize;

vx          = find(pos_list(:,1)>xmin & pos_list(:,1)<xmax & pos_list(:,2)>ymin & pos_list(:,2)<ymax);
ROIselect   = pos_list(vx,1:end);


% Plot the fiducials

close all
figure('Position',[100 200 400 400],'Name','Selected ROI')  
scatter(ROIselect(:,1),ROIselect(:,2),1);hold on; box on;
title('Selected ROI, Input for tracking');
xlabel('x [pxl]');
ylabel('y [pxl]');
axis([min(ROIselect(:,1)) max(ROIselect(:,1)) min(ROIselect(:,2)) max(ROIselect(:,2))])
set(gca,'FontSize',12)
    
%% 4. Track 
close all
% using the Crocker, Weeks, and Grier Algorithm (http://www.physics.emory.edu/rweeks/idl/index.html)

% mem   - number of time steps that a particle can be 'lost' and then recovered again
% dim   - default 2
% good  - eliminate if fewer than good valid positions
% quiet - 1 = no text

tic
param   = struct('mem',2,'dim',2,'good',1,'quiet',1);
res     = trackGT(ROIselect,1,param); % variable XYT, maximum displacement in unit of data

fprintf('\n -- Tracking done in %f sec--\n',toc)

% output res:
% 
% 1. X
% 2. Y
% 3. time
% 4. ID

%% 5. Filter and visualize tracks

close all, clc

min_length  = 5;
max_length  = 100; % max(peaks(:,frame))/10;
res_filt    = []; count = 0; track = [];

% Track length

for i=1:max(res(:,4))
     
    target = find(res(:,4)==i);
    
    if length(target)>min_length & length(target)<max_length;
    
    count       = count +1;
    track       = res(target,:);
    track(:,4)  = count;

    res_filt = vertcat(res_filt, track);
           
    else
    end
end

fprintf(['-- Tracks filtered, min length = ', num2str(min_length),' \n']);
fprintf(['-- Filtered tracks ', num2str(count),' - ',num2str(length(res_filt)/length(res)*100),' %% left \n']);

% Show dense regions and overlay with WF image

density = [];
NoNind  = [];

density(:,1)=res_filt(:,1);
density(:,2)=res_filt(:,2);

idx = rangesearch(density,density,1);
        
for i=1:length(idx);
NoNind(i,1)=length(idx{i,1});     % count the total number of neighbors for each point in dataset
end

figure('Position',[100 600 700 300],'Name',['Overlay'])
subplot(1,2,1);
imagesc(WF_image,[0 10*median(WF_image(:))]);hold on; axis square; title('STD + Points');
colormap(jet);
scatter(res_filt(:,1),res_filt(:,2),'r.');


subplot(1,2,2);
imagesc(WF_image,[0 10*median(WF_image(:))]);hold on;axis square;title('STD + Tracks');
colormap(gray);

for i=1:max(res_filt(:,4))
     
    target=find(res_filt(:,4)==i);
    
    if length(target)>min_length;
    
    plot(res_filt(target,1),res_filt(target,2));hold on;
    
    else
    end
end


% Track length histogram

tracklength = [];

for index = 1:max(res_filt(:,4))
    
    track       = find(res_filt(:,4)==index);
    tracklength = cat(1,tracklength,length(track));
    
    clear track
end

figure('Position',[100 100 750 300],'Name',['Track length'])
subplot(1,2,1)
hist(tracklength(tracklength>min_length),30);
title(['Mean Length = ',num2str(mean(tracklength))]);
xlabel('count');
ylabel('tracklength (frames)');
  
% Plot tracks longer than min_length

subplot(1,2,2)

for i=1:max(res_filt(:,4))
     
    target=find(res_filt(:,4)==i);
    
    if length(target)>min_length;
    
    plot(res_filt(target,1),res_filt(target,2));hold on;
    
    else
    end
end

    title(['All tracks longer than ',num2str(min_length) ' frames']);
    box on; axis square; axis equal;
    xlabel('pxl');
    ylabel('pxl');
   
%% 6. Export variable "tracks" to interact with @msdanalyzer

count = 1;       % counter to select field

tracks = [];

for i=1:max(res_filt(:,4))
     
    target = find(res_filt(:,4)==i);
        
    tracks{count,1}(:,1)=res_filt(target,3)*time_step;        % frame
    tracks{count,1}(:,2)=res_filt(target,1)*pxl_size/1e3;     % x in mum
    tracks{count,1}(:,3)=res_filt(target,2)*pxl_size/1e3;     % y in mum
    
    count=count+1;

end

cd(pathR);

filenameMSD=['Tracks_MSD_' name_base '.m'];
save(filenameMSD,'tracks');

fprintf('\n -- Done. Generated tracks for @msdanalyzer --\n')

% Export tracks as input for Hoze/inferenceMap

filenamec1=['Tracks_XYTID_' name_base '.m'];
filenamec2=['Tracks_Hoze_' name_base '.txt'];
filenamec3=['Tracks_InferenceMap_' name_base '.trxyt'];

cd(pathR);

% As Matlab without changes

save(filenamec1,'res_filt');

% 1 - x
% 2 - y
% 3 - time in seconds
% 4 - track ID

% For N. Hoze Script

forHoze=[];
forHoze(:,1)=res_filt(:,4);                 % index
forHoze(:,2)=res_filt(:,3)*time_step;       % time in sec
forHoze(:,3)=res_filt(:,1)*pxl_size/1e3;    % X in um
forHoze(:,4)=res_filt(:,2)*pxl_size/1e3;    % Y in um

dlmwrite(filenamec2,forHoze);

% For InferenceMap

forIM=[];
forIM(:,1)=res_filt(:,4);                   % ID
forIM(:,2)=res_filt(:,1)*pxl_size/1e3;      % X in um
forIM(:,3)=res_filt(:,2)*pxl_size/1e3;      % Y in um
forIM(:,4)=res_filt(:,3)*time_step;         % frame

dlmwrite(filenamec3,forIM,'delimiter','\t');

fprintf('\n -- Tracks Saved --\n')

