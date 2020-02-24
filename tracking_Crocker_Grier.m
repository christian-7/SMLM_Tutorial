%% sptPALM tracking and post-processing

% 1. Filter and vizualize
% 2. Select ROI
% 3. Track using Crocker and Grier
% 4. Filter tracks

clear, close all, clc, clear

%% %%%%%%%%%% Load Localizations %%%%%%%%%%

pathd = ('//Users/christian/Documents/Arbeit/SuperRes/20-02-21_M1_mEOS'); % Data path
pathR = ('//Users/christian/Documents/Arbeit/SuperRes/20-02-21_M1_mEOS'); % Results path

WF_image_name  = 'STD_A549_M1_mEos_30ms_gain300_c4.tif';
filename_peaks = 'A549_M1_mEos_30ms_gain300_c3_filt';     % filename of TS output file

pxl_size    = 160; % nm
time_step   = 0.03;

%%%%%%%%%% Load File %%%%%%%%%%

cd(pathd);
filename_peaks2=[filename_peaks '.csv'];
peaks=dlmread(filename_peaks2,',',1,0);
WF_image = imread(WF_image_name);

file = fopen(filename_peaks2);
line = fgetl(file);
h = regexp( line, ',', 'split' );

% Data from Matlab, Fangs software

x = strmatch('x [nm]',h);
y = strmatch('y [nm]',h);
frame = strmatch('frame',h);

% From Thunderstorm

x = strmatch('"x [nm]"',h);
y = strmatch('"y [nm]"',h);
frame = strmatch('"frame"',h);

fprintf('\n -- Data loaded --\n')

%% Filter and generate PosList

close all
min_int     = 3; % min photons
filtered    = []; density =[];
NNdist      = 200; % nm for the density plot
NoNind      = [];

%%%%% Filter

filter          = find(peaks(:,5)>min_int);
filtered        = peaks(filter,x);
filtered(:,2)   = peaks(filter,y);
filtered(:,3)   = peaks(filter,frame);

%%%%% Density

density=filtered(:,1:2);

idx = rangesearch(density,density,NNdist);
        
for i=1:length(idx);
NoNind(i,1)=length(idx{i,1});     % count the total number of neighbors for each point in dataset
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Position',[100 600 600 600])
imshow(flipud(WF_image));hold on;
colormap(hot)
scatter(density(:,1)/pxl_size,density(:,2)/pxl_size,1,NoNind);hold on;
box on;

figure('Position',[100 600 1200 400])
subplot(1,2,1)
scatter(filtered(:,1), filtered(:,2),1); hold on;
title([num2str(length(filtered)),' from ',num2str(length(peaks)),' left after filtering -> ', num2str(length(filtered)/length(peaks)*100),'  %',]);
axis([min(filtered(:,1)) max(filtered(:,1)) min(filtered(:,2)) max(filtered(:,2))])
xlabel('nm');
ylabel('nm');
box on; axis equal;

subplot(1,2,2)
scatter(density(:,1),density(:,2),5,NoNind,'o','filled');
colormap(jet)
title(['Localization density, NN Dist = ', num2str(NNdist) ' nm']);
axis([min(density(:,1)) max(density(:,1)) min(density(:,2)) max(density(:,2))])
xlabel('nm');
ylabel('nm');
whitebg(1,'k')
box on; axis equal;

% Create pos_list for track.m --> tracking in pxl/frames

pos_list=[];

pos_list(:,1)=filtered(:,1)/pxl_size;                   % in pxl
pos_list(:,2)=filtered(:,2)/pxl_size;                   % in pxl
pos_list(:,3)=filtered(:,3);                            % dt in frames

fprintf('\n -- Data filtered / PosList generated--\n')

%% Track using the Crocker, Weeks, and Grier Algorithm (http://www.physics.emory.edu/rweeks/idl/index.html)

% mem   - number of time steps that a particle can be 'lost' and then recovered again
% dim   - default 2
% good  - eliminate if fewer than good valid positions
% quiet - 1 = no text

tic
param=struct('mem',2,'dim',2,'good',1,'quiet',1);
res=trackGT(pos_list,3,param); % variable XYT, maximum displacement in unit of data

fprintf('\n -- Tracking done in %f sec--\n',toc)

% output res:
% 
% 1. X
% 2. Y
% 3. time
% 4. ID

%% Filter and visualize tracks

min_length  = 10;
max_length  = 100; %max(peaks(:,frame))/10;
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

figure('Position',[100 600 600 600])
imshow(flipud(WF_image));hold on;
colormap(jet)
scatter(density(:,1),density(:,2),1,NoNind);hold on;


% Track length histogram

tracklength = [];

for index = 1:max(res_filt(:,4))
    
    track       = find(res_filt(:,4)==index);
    tracklength = cat(1,tracklength,length(track));
    
    clear track
end

figure('Position',[100 500 600 200],'Name',['Track length'])
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

    title(['All tracks longer than ',num2str(min_length) 'frames']);
    xlabel('pxl');
    ylabel('pxl');
   
%% Export variable "tracks" to interact with @msdanalyzer

%%%%%%%%%%%%%%%%%%%%%% ROI selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%

upperx= max(res_filt(:,1));         
lowerx=0;

uppery= max(res_filt(:,2));
lowery= 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vx=find(res_filt(:,1) < upperx & res_filt(:,1) > lowerx & res_filt(:,2) < uppery & res_filt(:,2) > lowery);
res_wROI = res_filt(vx,:);

% figure
% scatter(res_wROI(:,1), res_wROI(:,2),2), hold on

fprintf('\n -- Done. Selected ROI --\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count       = 1;       % counter to select field

for i=1:max(res_wROI(:,4))
     
    target=find(res_wROI(:,4)==i);
    
    if length(target)>min_length;
    
    tracks{count,1}(:,1)=res_wROI(target,3)*time_step;     % frame
    tracks{count,1}(:,2)=res_wROI(target,1)*pxl_size/1e3;     % x in mum
    tracks{count,1}(:,3)=res_wROI(target,2)*pxl_size/1e3;     % y in mum
    
    count=count+1;
    
    else
    end
end

cd(pathR);

filenameMSD=['Tracks_MSD_' filename_peaks];
save(filenameMSD,'tracks');

fprintf('\n -- Done. Generated tracks for @msdanalyzer --\n')

% Export tracks as input for Hoze/inferenceMap

filenamec1=['Tracks_XYTID_' filename_peaks];
filenamec2=['Tracks_Hoze_' filename_peaks '.txt'];
filenamec3=['Tracks_InferenceMap_' filename_peaks '.trxyt'];

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

