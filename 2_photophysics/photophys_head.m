%% Analyze dye photophysics from SMLM data (tracking part)
% 
% Author: Christian Sieben, EPFL 
% sieben.christian@gmail.com
% April 2020

% 1. Load sample data
% 2. Track with Gap 0
% 3. Group/Merge the tracks (2D)
% 4. Track with Gap max(frames)

clear, clc, close all
SMLM_tutorial_main = '/Users/csi20/Documents/Arbeit/MatLab/SMLM_tutorial'; % for relative path

%% 1. Load sample data 

cd([SMLM_tutorial_main '/example_data/photophysics/A647']);

filename = 'A647_COT_1200mW_10ms_3_MMStack_1_Localizations_DC_Z';
locs = dlmread([filename '.csv'],',',1,0);

% Find the respective Columns

file   = fopen([filename '.csv']); % csv for TS
header = fgetl(file);
header = regexp(header, ',', 'split');

% ThunderStorm format

xCol                = strmatch('"x [nm]"',header);
yCol                = strmatch('"y [nm]"',header);
zCol                = strmatch('"z [nm]"',header);
photonsCol          = strmatch('"intensity [photons]"',header);
framesCol           = strmatch('"frame"',header);
LLCol               = strmatch('"loglikelihood"',header);

% Fangs format 

xCol                = strmatch('x [nm]',header);
yCol                = strmatch('y [nm]',header);
zCol                = strmatch('z [nm]',header);
photonsCol          = strmatch('intensity [photon]',header);
framesCol           = strmatch('frame',header);
LLCol               = strmatch('loglikelihood',header);
uncertaintyCol      = strmatch('uncertainty [nm]',header);

fprintf('\n -- Data loaded --\n')

figure
scatter(locs(:,xCol),locs(:,yCol),'.');
axis square, box on

%% 2. Track with Gap 0
%  track with a gap of 0 frames to combine blinking events spread over
%  several frames

max_disp = 50; % maximum displacement in unit of data
gap      = 0;  % number of time steps that a particle can be 'lost' and then recovered again
min_pos  = 1;  % eliminate if fewer than min_pos good valid positions
quiet    = 1;  % no text


pos_list(:,1) = locs(:,xCol);                   % in pxl
pos_list(:,2) = locs(:,yCol);                   % in pxl
pos_list(:,3) = locs(:,photonsCol);             % photons
pos_list(:,4) = locs(:,uncertaintyCol);         % uncertainty
pos_list(:,5) = locs(:,framesCol);              % dt in seconds


param = struct('mem',gap,'dim',2,'good',min_pos,'quiet',quiet);
res_gap0   = trackGT(pos_list,max_disp,param); % variable XYT, maximum displacement in pxl

save([filename '_tracks_Gap0'],'res_gap0');

fprintf('\n -- Tracking Done --\n')

%% 3. Group/Merge the tracks (2D)

% Initialize Variable for grouped locs

groupedx = []; groupedy = []; frame    = [];
groupedframe = []; groupedPhotons = []; groupedUncertainty=[];

% Photons = [];
% Uncertainty = [];
% subsetLL = [];

for index=1:max(res_gap0(:,end)); 
    
            vx = find(res_gap0(:,end)==index); % find the track ID
            
            if length(vx)<300 ; % select only tracks short tham 300 locs
            
            groupedx            = vertcat(groupedx,sum(res_gap0(res_gap0(:,6)==index,1))/length(res_gap0(res_gap0(:,6)==index,1)));
            groupedy            = vertcat(groupedy,sum(res_gap0(res_gap0(:,6)==index,2))/length(res_gap0(res_gap0(:,6)==index,2)));
            groupedframe        = vertcat(groupedframe, round(mean(res_gap0(res_gap0(:,6)==index,5))));
            groupedPhotons      = vertcat(groupedPhotons, sum(res_gap0(res_gap0(:,6)==index,3)));
            groupedUncertainty  = vertcat(groupedUncertainty, mean(res_gap0(res_gap0(:,6)==index,4)));
            
            else end % else its a bead
end

locs_grouped = [];
locs_grouped = [groupedx, groupedy, groupedPhotons, groupedUncertainty, groupedframe];

% save([filename '_tracks_Gap0_Merged'],'locs_grouped');

fprintf('\n -- 3. Tracks Merged and Saved --\n')

%% 3. Tracking with Gap max(frames)
%  track with a gap of max frames to combine blinking events into molecules

% Set the grouping parameters

max_disp = 50;                      % maximum displacement in unit of data
gap      = max(locs_grouped(:,5));  % number of time steps that a particle can be 'lost' and then recovered again
min_pos  = 1;                       % eliminate if fewer than min_pos good valid positions
quiet    = 1;                       % no text

pos_list = [];

pos_list(:,1) = locs_grouped(:,1); % x
pos_list(:,2) = locs_grouped(:,2); % y
pos_list(:,3) = locs_grouped(:,3); % Photons
pos_list(:,4) = locs_grouped(:,4); % Uncertainty
pos_list(:,5) = locs_grouped(:,5); % Frame

pos_list = sortrows(pos_list,5);

param   = struct('mem',gap,'dim',2,'good',min_pos,'quiet',quiet);
res_gapMax     = trackGT(pos_list,max_disp,param); % variable XYT, maximum displacement in pxl

fprintf('\n -- Tracking Done --\n')

save([filename '_tracks_GapMax' ],'res_gapMax');

fprintf('\n -- 3. Tracks Saved --\n')

