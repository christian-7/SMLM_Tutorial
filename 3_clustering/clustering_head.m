%% Perform clustering of SMLM data
% 
% Author: Christian Sieben, EPFL 
% sieben.christian@gmail.com
% April 2020

% 1. Load sample data
% 2. Calculate Ripley's K statistics
% 3. Perform Ripley's K clustering

% Ripleys K
% This approach compares the measured distribution of single-molecule localizations to a simulated random distribution, and provides information whether clustering occurs. 
% The maximum of the H function reflects the average size of clusters. 
% The amplitude of the H-function is a measure of the degree of clustering. 
% The variance of the amplitude reflects the variance of the spatial organization of proteins within a cluster.


clear, clc, close all
SMLM_tutorial_main = '/Users/christian/Documents/Arbeit/MatLab/SMLM_tutorial';

%% 1. Load sample data

% cd([SMLM_tutorial_main '/example_data/clustering']);
% filename = 'A549_EGFR_Co_800ms_10ms_3_MMStack_locResults.dat';

folder = '/Volumes/Seagate Expansion Drive/Christian_lebpc4/Inflammasome_backup/2020-03-03_CS_Inflammasome/locResultsTresh100/Sample1_Unstimulated_18_1';
cd(folder)
filename = 'Sample1_Unstimulated_18_1_MMStack_1_Localizations';

locs=dlmread([filename '.csv'],',',1,0);

% Find the respective Columns

file = fopen([filename '.csv']); % csv for TS
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


fprintf('\n -- Data loaded --\n')

%% 2. Peform DBSCAN clustering (with manual ROI)

cd([SMLM_tutorial_main '/3_clustering']);

pxlsize = 200;

k       = 5; 
Eps     = 50;

[subset, FiC, CpA, LpC] = DBSCAN_with_ROI(locs, xCol, yCol, pxlsize, k, Eps);

%% 3. Calculate Ripley's K statistics (with manual ROI)
% Generate loc file for LAMA *_Lama.txt

close all

[locs_ROI] = manualROI(locs,xCol,yCol,200);

ROI     = [min(locs_ROI(:,xCol)) max(locs_ROI(:,xCol)) ... 
           min(locs_ROI(:,yCol)) max(locs_ROI(:,yCol))];

maxk    = 1000;    % maximum distance
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

% cd(folder)
% saveAsMALK(locs_ROI,xCol, yCol, framesCol, photonsCol,filename)

%% 4. Calculate Ripley's K and L statistics (with manual ROI)

% Code from: Shivanandan A, et al (2015) PLOS ONE 10(3): e0118767. https://doi.org/10.1371/journal.pone.0118767
% https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0118767

close all

radius = 5:5:500;

[locs_ROI] = manualROI(locs,xCol,yCol,200);

ROI     = [min(locs_ROI(:,xCol)) max(locs_ROI(:,xCol)) ... 
           min(locs_ROI(:,yCol)) max(locs_ROI(:,yCol))];

ROI_cropped = find(locs_ROI(:,xCol)>ROI(1)+max(radius) & locs_ROI(:,xCol)<ROI(2)-max(radius) & ...
                   locs_ROI(:,yCol)>ROI(3)+max(radius) & locs_ROI(:,yCol)<ROI(4)-max(radius));


[TrueLrr, TrueK] = calculateLKfunction(locs_ROI(:,1:2), radius);

[max_value,index] = max(TrueLrr)

figure('Position',[100 400 600 300])

subplot(1,2,1)
plot(radius,TrueK,'-b');hold on
ylabel('K(r)'); box on; axis square; xlabel('distance [nm]');

subplot(1,2,2)
plot(radius,TrueLrr,'b');hold on
title(['Maximum at ' num2str(max(radius(index))) ' nm'])
ylabel('L(r)-r');box on; axis square; xlabel('distance [nm]');

% Var1(:,1) = radius;
% Var1(:,3) = TrueLrr;

% cd(folder)
% saveAsMALK(locs_ROI,xCol, yCol, framesCol, photonsCol, [filename '_ROI2_'])

%% 5. Local linearised Ripley's K density estimator
% Code provided by Juliette Griffie, EPFL (see Lr.m)

[locs_ROI] = manualROI(locs,xCol,yCol,200);

ROI     = [min(locs_ROI(:,xCol)) max(locs_ROI(:,xCol)) ... 
           min(locs_ROI(:,yCol)) max(locs_ROI(:,yCol))];
       
radius  = 100;

ROI_cropped = find(locs_ROI(:,xCol)>ROI(1)+radius & locs_ROI(:,xCol)<ROI(2)-radius & ...
                   locs_ROI(:,yCol)>ROI(3)+radius & locs_ROI(:,yCol)<ROI(4)-radius);
               
locs_ROI_cropped = locs_ROI(ROI_cropped, xCol:yCol);

figure('Position',[100 400 300 300])
scatter(locs_ROI(:,xCol),locs_ROI(:,yCol),'k.');hold on;
scatter(locs_ROI_cropped(:,xCol),locs_ROI_cropped(:,yCol),'r.');hold on;
axis square; box on; xlabel('x (nm)'); ylabel('y (nm)');


D = pdist2(locs_ROI_cropped,locs_ROI(:,xCol:yCol));
L = [];
Area = (ROI(2)-ROI(1))*(ROI(4)-ROI(3));

for i = 1:size(locs_ROI_cropped,1);
    delta = D(i,:);
    delta(delta>radius)=[];
%     L = [L; sqrt((Area/(pi*size(locs_ROI_cropped,1)))*(size(delta,2)-1))];
    L=[L; sqrt((Area/(pi*size(locs_ROI_cropped,1)))*(size(delta,2)-1))-radius];
    
end

figure('Position',[400 400 300 300])
scatter(locs_ROI_cropped(:,1),locs_ROI_cropped(:,2),5, L,'filled')
axis square;box on; hold on; xlabel('x (nm)'); ylabel('y (nm)');

%% 6. Calculate L Curve (same ROI as 4)
% Code provided by Juliette Griffie, EPFL (see Lr.m)

radius = 5:5:500;

[locs_ROI] = manualROI(locs,xCol,yCol,200);

ROI_cropped = find(locs_ROI(:,xCol)>ROI(1)+max(radius) & locs_ROI(:,xCol)<ROI(2)-max(radius) & ...
                   locs_ROI(:,yCol)>ROI(3)+max(radius) & locs_ROI(:,yCol)<ROI(4)-max(radius));
               
locs_ROI_cropped = locs_ROI(ROI_cropped, xCol:yCol);

D           = pdist2(locs_ROI_cropped,locs_ROI(:,xCol:yCol));
D_vector    = D(:);
D_vector(D_vector==0)=[];           % remove zeros
D_vector(D_vector>max(radius))=[];  % remove all values larger than max(radius)

figure('Position',[100 400 300 300])
H = histogram(D_vector, size(radius,2));
box on
axis square
xlabel('distance (nm)'); ylabel('counts')

count = H.Values;

L_curves = [];
L_curves = [L_curves,0];

for i=1:size(radius,2)
    
    L_curves=[L_curves, sqrt((Area/(pi*(size(locs_ROI_cropped,1)^2)))*sum(count(1:i)))-radius(i)];
    
end

radius = [0,radius];
[max_value,index]=max(TrueLrr);

figure('Position',[400 400 300 300])
plot(radius,L_curves);
box on
axis square
xlabel('radius (nm)'); ylabel('L(r)-r')
title(['Maximum at ' num2str(max(radius(index))) ' nm'])
