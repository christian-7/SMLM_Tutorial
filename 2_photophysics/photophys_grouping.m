%% Grouping/merging by using 
% 
% Author: Christian Sieben, HZI 
% sieben.christian@gmail.com
% Sep 2020

% 1.    Load Locs
% 2.    Define grouping parameters, perfomr grouping
% 3. 	Merge and save localizations

clear, clc, close all

SMLM_tutorial_main = '/Users/csi20/Documents/Arbeit/MatLab/SMLM_tutorial'; % for relative path

%% 1. Load locs

cd([SMLM_tutorial_main '/example_data/photophysics/A647']);

filename = 'A647_COT_1200mW_10ms_3_MMStack_1_Localizations_DC_Z';
locs = dlmread([filename '.csv'],',',1,0);

% Find the respective columns

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


fprintf('\n -- 1. Data loaded --\n')


%% 2. Perform grouping 

% Set the grouping parameters

max_disp = 30;    % maximum displacement, e.g.2*sigma
gap      = 1344;  % number of time steps that a particle can be 'lost' and then recovered again, e.g.  mean dark time
min_pos  = 1;     % eliminate if fewer than min_pos good valid positions
quiet    = 1;     % no text

pos_list(:,1) = locs(:,xCol);                   % in nm
pos_list(:,2) = locs(:,yCol);                   % in nm
pos_list(:,3) = locs(:,photonsCol);             % photons
pos_list(:,4) = locs(:,uncertaintyCol);         % uncertainty
pos_list(:,5) = locs(:,framesCol);              % frames

param = struct('mem',gap,'dim',2,'good',min_pos,'quiet',quiet);
res   = trackGT(pos_list,max_disp,param); % variable XYT, maximum displacement in pxl

fprintf('\n -- 2. Tracking Done --\n')

%% 3. Merge localizations

% Initialize Variables for grouped locs

groupedx = []; groupedy = []; frame    = [];
groupedframe = []; groupedPhotons = []; groupedUncertainty=[];

for index = 1:max(res(:,end)); % for all tracks 
    
            vx = find(res(:,end)==index); % find the track ID
            
            if length(vx)<300 ; % select only tracks short tham 300 locs, to filter Au fiducials
            
            groupedx            = vertcat(groupedx,sum(res(res(:,6)==index,1))/length(res(res(:,6)==index,1)));
            groupedy            = vertcat(groupedy,sum(res(res(:,6)==index,2))/length(res(res(:,6)==index,2)));
            groupedframe        = vertcat(groupedframe, round(mean(res(res(:,6)==index,5))));
            groupedPhotons      = vertcat(groupedPhotons, sum(res(res(:,6)==index,3)));
            groupedUncertainty  = vertcat(groupedUncertainty, mean(res(res(:,6)==index,4)));
            
            else end % else its a bead
end

locs_grouped = [];
locs_grouped = [groupedx, groupedy, groupedPhotons, groupedUncertainty, groupedframe];

% Plot figure to visualize before/after merging

figure
scatter(locs(:,xCol), locs(:,yCol),'b.'); hold on;
scatter(locs_grouped(:,1), locs_grouped(:,2),15,'r')

% Save merged localizations

% save(['Locs_Merged_' filename],'locs_grouped');

fprintf('\n -- 3. Locs Merged and Saved --\n')
