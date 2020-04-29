%% Load HTP Dataset and perform tracking 

% unsing the Crocker, Weeks, and Grier Algorithm (http://www.physics.emory.edu/rweeks/idl/index.html)

% clear, close all, clc, clear
%%  %%%%%%%%%% INPUT Parameters %%%%%%%%%%
function [res]=CG_tracking_TS_input(filename_peaks, max_disp, gap, saveYN);

% filename_peaks='NB_Nbr5_A647_COT_1500mW_1_MMStack_locResults_processed';     % filename of TS output file
% 
% max_disp = 20;   % in pxl
min_pos = 1;     % good - eliminate if fewer than good valid positions
% gap = 0;         % mem - number of time steps that a particle can be 'lost' and then recovered again
quiet = 1;       % quiet - 1 = no text

%% Load Data and find the columns

peaks=dlmread([filename_peaks '.csv'],',',1,0);

file = fopen([filename_peaks '.csv']);
line = fgetl(file);
h = regexp( line, ',', 'split' );

x           = strmatch('"x [nm]"',h);
y           = strmatch('"y [nm]"',h);
frame       = strmatch('"frame"',h);
photons     = strmatch('"intensity [photon]"',h);
sigma       = strmatch('"sigma [nm]"',h);
uncertainty = strmatch('"uncertainty_xy [nm]"',h);

% Create pos_list for track.m

pos_list(:,1) = peaks(:,x);                   % in pxl
pos_list(:,2) = peaks(:,y);                   % in pxl
pos_list(:,3) = peaks(:,photons);             % photons
pos_list(:,4) = peaks(:,uncertainty);         % uncertainty
pos_list(:,5) = peaks(:,frame);               % dt in seconds

fprintf('\n -- Data Loaded --\n')

%% Track unsing the Crocker, Weeks, and Grier Algorithm (http://www.physics.emory.edu/rweeks/idl/index.html)

param = struct('mem',gap,'dim',2,'good',min_pos,'quiet',quiet);
res   = trackGT(pos_list,max_disp,param); % variable XYT, maximum displacement in pxl

fprintf('\n -- Tracking Done --\n')

% res1 = x
% res2 = y
% res3 = photons
% res4 = uncertainty
% res5 = frame
% res6 = track ID

%% Plot All Tracks

% scatter(res(:,1),res(:,2),5,res(:,6));
% 
% for i=1:max(res(:,6))
%     
%     vx=find(res(:,6)==i);
%     
%     plot(res(vx,1),res(vx,2));hold on;
% 
% end

%% Save tracks

if saveYN==1;

filenamec1=['Tracks_MaxGap' filename_peaks];

save(filenamec1,'res');

fprintf('\n -- Tracks Saved --\n')

else
    
fprintf('\n -- Tracks Not Saved --\n')
    
end

% 1 - x
% 2 - y
% 3 - time in seconds
% 4 - track ID
% 5 - time in seconds
% 6 - track ID

end