%% sptPALM tracking and post-processing

% Author: Christian Sieben, EPFL 
% sieben.christian@gmail.com
% January 2020

% Content

% 1. Load loc file 
% 2. Filter and vizualize
% 3. Select ROI
% 4. Track
% 5. Filter tracks
% 6. Export tracks

%% 1. Load Localizations

clear, close all, clc, clear

pathd = ('//Users/christian/Documents/Arbeit/SuperRes/20-02-21_M1_mEOS'); % Data path
pathR = ('//Users/christian/Documents/Arbeit/SuperRes/20-02-21_M1_mEOS'); % Results path

WF_image_name  = 'STD_A549_M1_mEos_30ms_gain300_c4.tif';
filename_peaks = 'A549_M1_mEos_30ms_gain300_c4_filt';     % filename of TS output file

pxl_size    = 160;  % nm
time_step   = 0.03; % in sec

%%%%%%%%%% Load File %%%%%%%%%%

cd(pathd);
filename_peaks2=[filename_peaks '.csv'];
peaks=dlmread(filename_peaks2,',',1,0);
WF_image = imread(WF_image_name);

file = fopen(filename_peaks2);
line = fgetl(file);
h = regexp( line, ',', 'split' );

% Data from Matlab (Fangs, fitcspline)

x = strmatch('x [nm]',h);
y = strmatch('y [nm]',h);
frame = strmatch('frame',h);

% From Thunderstorm

x = strmatch('"x [nm]"',h);
y = strmatch('"y [nm]"',h);
frame = strmatch('"frame"',h);

fprintf('\n -- Data loaded --\n')

%% 2. Filter and generate PosList

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

%% 3. Select ROI
close all
pxlsize = 1;

heigth  = round((max(pos_list(:,2))-min(pos_list(:,2)))/pxlsize);
width   = round((max(pos_list(:,1))-min(pos_list(:,1)))/pxlsize);
im      = hist3([pos_list(:,1),pos_list(:,2)],[width heigth]); % heigth x width

% Select rectangles

rect = [];

figure('Position',[100 200 400 400])
imagesc(imrotate(im,90));
title('Please select ROI');
xlabel('x [pxl]');
ylabel('y [pxl]');

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
scatter(ROIselect(:,1),ROIselect(:,2),1);hold on;
title('Selected ROI, Input for tracking');
xlabel('x [pxl]');
ylabel('y [pxl]');
    
%% 4. Track 

% using the Crocker, Weeks, and Grier Algorithm (http://www.physics.emory.edu/rweeks/idl/index.html)

% mem   - number of time steps that a particle can be 'lost' and then recovered again
% dim   - default 2
% good  - eliminate if fewer than good valid positions
% quiet - 1 = no text

tic
param=struct('mem',2,'dim',2,'good',1,'quiet',1);
res=trackGT(ROIselect,3,param); % variable XYT, maximum displacement in unit of data

fprintf('\n -- Tracking done in %f sec--\n',toc)

% output res:
% 
% 1. X
% 2. Y
% 3. time
% 4. ID

%% 5. Filter and visualize tracks

close all, clc

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

figure('Position',[500 600 600 200],'Name',['Track length'])
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

