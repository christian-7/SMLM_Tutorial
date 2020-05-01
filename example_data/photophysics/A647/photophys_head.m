%% Analyze dye photophysics from SMLM data
% 
% Author: Christian Sieben, EPFL 
% sieben.christian@gmail.com
% April 2020

% 1. Load sample data
% 2. Track with Gap 0
% 3. Group/Merge the tracks (2D)
% 4. Track with Gap max(frames)

clear, clc, close all
SMLM_tutorial_main = '/Users/christian/Documents/Arbeit/MatLab/SMLM_tutorial';

%% 1. Load sample data 

cd([SMLM_tutorial_main '/example_data/photophysics/A647']);

% filename = 'HAmOrange_NB_41_1_locs_ROI.csv';
filename = 'A647_COT_1200mW_10ms_3_MMStack_1_Localizations_DC_Z.csv';
locs=dlmread(filename,',',1,0);

% Find the respective Columns

file = fopen(filename); % csv for TS
header = fgetl(file);
header = regexp(header, ',', 'split');

xCol                = strmatch('"x [nm]"',header);
yCol                = strmatch('"y [nm]"',header);
zCol                = strmatch('"z [nm]"',header);
photonsCol          = strmatch('"intensity [photons]"',header);
framesCol           = strmatch('"frame"',header);
LLCol               = strmatch('"loglikelihood"',header);

xCol                = strmatch('x [nm]',header);
yCol                = strmatch('y [nm]',header);
zCol                = strmatch('z [nm]',header);
photonsCol          = strmatch('intensity [photon]',header);
framesCol           = strmatch('frame',header);
LLCol               = strmatch('loglikelihood',header);
uncertaintyCol      = strmatch('uncertainty [nm]',header);


fprintf('\n -- Data loaded --\n')

%% 2. Track with Gap 0
%  track with gap 0 to combine blinking events

max_disp = 50;
gap      = 0;
min_pos  = 1;
quiet    = 1;


pos_list(:,1) = locs(:,xCol);                   % in pxl
pos_list(:,2) = locs(:,yCol);                   % in pxl
pos_list(:,3) = locs(:,photonsCol);             % photons
pos_list(:,4) = locs(:,uncertaintyCol);         % uncertainty
pos_list(:,5) = locs(:,framesCol);               % dt in seconds


param = struct('mem',gap,'dim',2,'good',min_pos,'quiet',quiet);
res   = trackGT(pos_list,max_disp,param); % variable XYT, maximum displacement in pxl

fprintf('\n -- Tracking Done --\n')

%% 3. Group/Merge the tracks (2D)

groupedx=[];
groupedy=[];
frame=[];
groupedframe=[];
groupedPhotons=[];
Photons=[];
Uncertainty = [];
subsetLL=[];
groupedUncertainty=[];


for index=1:max(res(:,6)); 
            
            if length(vx)<300 ; % select only tracks short tham 300 locs
            
            groupedx            = vertcat(groupedx,sum(res(res(:,6)==index,1))/length(res(res(:,6)==index,1)));
            groupedy            = vertcat(groupedy,sum(res(res(:,6)==index,2))/length(res(res(:,6)==index,2)));
            groupedframe        = vertcat(groupedframe, round(mean(res(res(:,6)==index,5))));
            groupedPhotons      = vertcat(groupedPhotons, sum(res(res(:,6)==index,3)));
            groupedUncertainty  = vertcat(groupedUncertainty, mean(res(res(:,6)==index,4)));
            
            else end % else its a bead
end

locs_grouped = [];
locs_grouped = [groupedx, groupedy, groupedPhotons, groupedUncertainty, groupedframe];


filenamec2=['Tracks_Gap0_Merged_' filename];
% save(filenamec2,'subsetLL');

fprintf('\n -- 3. Tracks Merged and Saved --\n')

%% 3. Tracking with Gap max(frames)

% Set the grouping parameters

max_disp    = 50;
gap         = max(locs_grouped(:,5));
min_pos     = 1;
quiet       = 1;

pos_list = [];

pos_list(:,1) = locs_grouped(:,1); % x
pos_list(:,2) = locs_grouped(:,2); % y
pos_list(:,3) = locs_grouped(:,3); % Photons
pos_list(:,4) = locs_grouped(:,4); % Uncertainty
pos_list(:,5) = locs_grouped(:,5); % Frame

pos_list = sortrows(pos_list,5);

param   = struct('mem',gap,'dim',2,'good',min_pos,'quiet',quiet);
res     = trackGT(pos_list,max_disp,param); % variable XYT, maximum displacement in pxl

fprintf('\n -- Tracking Done --\n')

filenamec3=['Tracks_GapMax_' filename];

% save(filenamec3,'res');

fprintf('\n -- 3. Tracks Saved --\n')



