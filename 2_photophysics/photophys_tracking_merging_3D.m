clear, clc, close all

%% Load the localization file

filename='A647_COT_800mW_10ms_1_MMStack_1_Localizations_DC_Z';

% 1. Tracking with Gap 0

[res]=CG_tracking(filename,50,0,0); % track with gap 0

%[res]=CG_tracking_TS_input(filename,50,0,0); % track with gap 0 to combine blinking events

filenamec1=['Tracks_Gap0_' filename];

% save(filenamec1,'res');

fprintf('\n -- 1. Tracks Saved --\n')

% 2. Merge 

% res1 = x
% res2 = y
% res3 = photons
% res4 = uncertainty
% res5 = frame
% res6 = track ID

groupedx=[];
groupedy=[];
groupedz=[];
frame=[];
groupedframe=[];
groupedID=[];
groupedPhotons=[];
Photons=[];
Uncertainty = [];
subsetLL=[];
groupedUncertainty=[];


for index=1:max(res(:,6)); 
    
            vx=find(res(:,6)==index);
            
            if length(vx)<300 ; % select only tracks short tham 300 locs
    
            clusterx = [];
            clustery = [];
            clusterz = [];
            clusterxC=[];
            clusteryC=[];
            clusterzC=[];
                                                   
            clusterx = res(vx,1);
            clustery = res(vx,2);
            clusterz = res(vx,3);
            frame = res(vx,5);
            
            clusterxC=sum(clusterx)/length(clusterx);
            clusteryC=sum(clustery)/length(clustery);
            clusterzC=sum(clusterz)/length(clusterz);
            Photons=sum(res(vx,3));
            Uncertainty = mean(res(vx,4));
            
            
            groupedx            = vertcat(groupedx,clusterxC);
            groupedy            = vertcat(groupedy,clusteryC);
            groupedz            = vertcat(groupedz,clusterzC);
            groupedframe        = vertcat(groupedframe, round(mean(frame)));
            groupedID           = vertcat(groupedID, index);
            groupedPhotons      = vertcat(groupedPhotons, Photons);
            groupedUncertainty  = vertcat(groupedUncertainty, Uncertainty);
            
            
            else end % else its a bead
end


subsetLL(:,1) = groupedx;
subsetLL(:,2) = groupedy;
subsetLL(:,3) = groupedz;
subsetLL(:,4) = groupedPhotons;
subsetLL(:,5) = groupedUncertainty;
subsetLL(:,6) = groupedframe;
subsetLL(:,7) = groupedID;

filenamec2=['Tracks_Gap0_Merged_' filename];
% save(filenamec2,'subsetLL');

fprintf('\n -- 2. Tracks Merged and Saved --\n')

%  3. Tracking with Gap max(frames)

% Set the grouping parameters

gap         = max(subsetLL(:,4));
min_pos     = 1;
max_disp    = 50;
quiet       = 1;

pos_list=[];

pos_list(:,1)=subsetLL(:,1); % x
pos_list(:,2)=subsetLL(:,2); % y
pos_list(:,3)=subsetLL(:,3); % z
pos_list(:,4)=subsetLL(:,4); % Photons
pos_list(:,5)=subsetLL(:,5); % Uncertainty
pos_list(:,6)=subsetLL(:,6); % Frame

pos_list=sortrows(pos_list,6);

param=struct('mem',gap,'dim',2,'good',min_pos,'quiet',quiet);
res = trackGT(pos_list,max_disp,param); % variable XYT, maximum displacement in pxl

fprintf('\n -- Tracking Done --\n')

filenamec3=['Tracks_GapMax_' filename];

save(filenamec3,'res');

fprintf('\n -- 3. Tracks Saved --\n')

clear

%% Test the merging parameters

% Load the dataset again

filename='2016-08-19_NB_A647_PLL-PM-glas_COT_10ms_1000mW_4_MMStack_Pos0_locResults_DC';

% Track with the new merging parameters

[res]=CG_tracking(filename,40,8000,0); % track with gap time (i.e. 3000 frames) and group radius (i.e. 50 nm) 

% Merge into molecule positions

% res1 = x
% res2 = y
% res3 = photons
% res4 = uncertainty
% res5 = frame
% res6 = track ID

groupedx=[];
groupedy=[];
frame=[];
groupedframe=[];
groupedID=[];
groupedPhotons=[];
Photons=[];
Uncertainty = [];
subsetLL=[];
groupedUncertainty=[];


for index=1:max(res(:,6)); 
    
            vx=find(res(:,6)==index);
            
            if length(vx)<1000 ; % select only tracks short tham 300 locs
    
            clusterx=[];
            clustery=[];
            clusterxC=[];
            clusteryC=[];
                                                   
            clusterx=res(vx,1);
            clustery=res(vx,2);
            frame=res(vx,5);
            
            clusterxC=sum(clusterx)/length(clusterx);
            clusteryC=sum(clustery)/length(clustery);
            Photons=sum(res(vx,3));
            Uncertainty = mean(res(vx,4));
            
            
            groupedx            = vertcat(groupedx,clusterxC);
            groupedy            = vertcat(groupedy,clusteryC);
            groupedframe        = vertcat(groupedframe, round(mean(frame)));
            groupedID           = vertcat(groupedID, index);
            groupedPhotons      = vertcat(groupedPhotons, Photons);
            groupedUncertainty  = vertcat(groupedUncertainty, Uncertainty);
            
            
            else end % else its a bead
end


subsetLL(:,1) = groupedx;
subsetLL(:,2) = groupedy;
subsetLL(:,3) = groupedPhotons;
subsetLL(:,4) = groupedUncertainty;
subsetLL(:,5) = groupedframe;
subsetLL(:,6) = groupedID;

% 3. Tracking with Gap max(frames)

% Set the grouping parameters

gap         = max(subsetLL(:,4));
min_pos     = 1;
max_disp    = 40;
quiet       = 1;

pos_list=[];

pos_list(:,1)=subsetLL(:,1); % x
pos_list(:,2)=subsetLL(:,2); % y
pos_list(:,3)=subsetLL(:,3); % Photons
pos_list(:,4)=subsetLL(:,4); % Uncertainty
pos_list(:,5)=subsetLL(:,5); % Frame

pos_list=sortrows(pos_list,5);

param=struct('mem',gap,'dim',2,'good',min_pos,'quiet',quiet);
res = trackGT(pos_list,max_disp,param); % variable XYT, maximum displacement in pxl

%  go to photophys_processing script
