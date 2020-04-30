%% Perform clustering of SMLM data
% 
% Author: Christian Sieben, EPFL 
% sieben.christian@gmail.com
% April 2020

% 1. Load sample data
% 2. Calculate Ripley's K statistics
% 3. Perform Ripley's K clustering

clear, clc, close all
SMLM_tutorial_main = '/Users/christian/Documents/Arbeit/MatLab/SMLM_tutorial';

%% 1. Load sample data

cd([SMLM_tutorial_main '/example_data/clustering']);

% filename = 'HAmOrange_NB_41_1_locs_ROI.csv';
filename = 'A549_EGFR_Co_800ms_10ms_3_MMStack_locResults.dat';
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
photonsCol          = strmatch('intensity [photons]',header);
framesCol           = strmatch('frame',header);
LLCol               = strmatch('loglikelihood',header);


fprintf('\n -- Data loaded --\n')

%% 2. Peform clustering

cd([SMLM_tutorial_main '/3_clustering']);

pxlsize = 200;

k       = 5; 
Eps     = 50;

[subset, FiC, CpA, LpC] = DBSCAN_with_ROI(locs, xCol, yCol, pxlsize, k, Eps);

%% 3. Calculate Ripley's K statistics

[locs_ROI] = manualROI(locs,xCol,yCol,200);

ROI     = [min(locs_ROI(:,xCol)) max(locs_ROI(:,xCol)) ... 
           min(locs_ROI(:,yCol)) max(locs_ROI(:,yCol))];

maxk    = 300;    % maximum distance
stepk   = 10;     % step size

[K, K_rand] = calculateRipleysK(locs_ROI, xCol, yCol, ROI, maxk, stepk);

L = sqrt(K(:,2)/pi);
H = L - K(:,1);
M = K(:,2)./(pi*(K(:,1).^2));

L_rand = sqrt(K_rand(:,2)/pi);
H_rand = L_rand - K_rand(:,1);
M_rand = K_rand(:,2)./(pi*(K_rand(:,1).^2));

figure('Position',[100 400 500 500])

subplot(2,2,1)
plot(K(:,1),K(:,2),'-b');hold on
plot(K_rand(:,1),K_rand(:,2),'-r');hold on
ylabel('K(r)'); box on; axis square; xlabel('distance [nm]');

subplot(2,2,2)
plot(K(:,1),L,'b');hold on
plot(K_rand(:,1),L_rand,'r');hold on
ylabel('L(r)');box on; axis square; xlabel('distance [nm]');

subplot(2,2,3)
plot(K(:,1),H,'b');hold on
plot(K_rand(:,1),H_rand,'r');hold on
ylabel('L(r)-r');box on; axis square; xlabel('distance [nm]');

subplot(2,2,4)
plot(K(:,1),M,'b');hold on
plot(K_rand(:,1),M_rand,'r');hold on
ylabel('K(r)/(\pi*r^2)'); box on; axis square; xlabel('distance [nm]');

%% 4. Localization dependent Ripley's K
