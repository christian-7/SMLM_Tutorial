%% Load the data

clear, clc, close all

locname = 'FOV_DL755_10ms_500mW_7_MMStack_1_Localizations_DC_Z';

locs        = dlmread([locname '.csv'],',',1,0);

file = fopen([locname '.csv']);
line = fgetl(file);
header1 = regexp( line, ',', 'split' );

xCol            = strmatch('x [nm]',header1);
yCol            = strmatch('y [nm]',header1);
zCol            = strmatch('z [nm]',header1);
framesCol       = strmatch('frame',header1);
LLCol           = strmatch('loglikelihood',header1);
photonsCol      = strmatch('intensity [photon]',header1);
uncertaintyCol  = strmatch('uncertainty [nm]',header1);

%% DBSCAN the data

dataDBS      = [];
dataDBS(:,1) = locs(:,xCol); % x in mum
dataDBS(:,2) = locs(:,yCol); % y in mum

% Run DBSCAN 

k   = 5;                                                    % minimum number of neighbors within Eps
Eps = 80;                                                    % minimum distance between points, nm

[class,type]=DBSCAN(dataDBS,k,Eps);                         % uses parameters specified at input
class2=transpose(class);                                    % class - vector specifying assignment of the i-th object to certain cluster (m,1)
type2=transpose(type);                                      % (core: 1, border: 0, outlier: -1)

coreBorder = [];
coreBorder = find(type2 >= 0);

subset          = [];
subset          = locs(coreBorder,1:end);
subset(:,end+1) = class2(coreBorder);

subsetP = [];
subsetP(:,1)    = dataDBS(coreBorder,1);
subsetP(:,2)    = dataDBS(coreBorder,2);
subsetP(:,3)    = class2(coreBorder);

identifiedClusters = max(class2(coreBorder))

%% Save the identified clusters

cd('Z:\Christian-Sieben\data_HTP\2017-07-24_TCI_Calibration_Alice');

save([locname '_DBSCAN_res.mat'],'subset');

%% Plot DBSCAN results 

figure('Position',[100 600 600 300]) % all data from both channels
subplot(1,2,1)
scatter(dataDBS(:,1),dataDBS(:,2),1);
axis([min(dataDBS(:,1)) max(dataDBS(:,1)) min(dataDBS(:,2)) max(dataDBS(:,2))])
box on

subplot(1,2,2)
scatter(subsetP(:,1),subsetP(:,2),1,mod(subsetP(:,3),10))
title('identified Clusters')
axis on
axis([min(dataDBS(:,1)) max(dataDBS(:,1)) min(dataDBS(:,2)) max(dataDBS(:,2))])
box on

%% Overlay the molecules

i = 4;

target = find(subset(:,end)==i);

selectedMol = subset(target,:);

figure
scatter3(selectedMol(:,xCol),selectedMol(:,yCol),selectedMol(:,zCol));
title(num2str(length(selectedMol)));
axis equal


%% Calculate Track Spread

%%%%%%%%%%%%
% The lateral spread of localizations 
%%%%%%%%%%%%


allclustersCx=[];
allclustersCy=[];
allclustersCz=[];
                                               
for index=1:max(subset(:,end)); % for all tracks
    
    vx = find(subset(:,end)==index);
    
    clusterx = [];
    clustery = [];
    clusterz = []; 
    
        if  length(vx)>1;                                                             % if nan, copy  frame number in new

            clusterx=subset(vx,xCol);
            clustery=subset(vx,yCol);
            clusterz=subset(vx,zCol);
            
            clusterxC=sum(clusterx)/length(clusterx);
            clusteryC=sum(clustery)/length(clustery);
            clusterzC=sum(clusterz)/length(clusterz);
            clusterx=clusterx-clusterxC;
            clustery=clustery-clusteryC;
            clusterz=clusterz-clusterzC;
            
            if clusterx>50; 
            else
            

            allclustersCx=vertcat(allclustersCx,clusterx);
            allclustersCy=vertcat(allclustersCy,clustery);
            allclustersCz=vertcat(allclustersCz,clusterz);
            
            end
        else end

end


figure('Position',[100 100 500 500],'name','Spread of Molecule (Overlay)')

subplot(3,2,[1 3 5])
scatter3(allclustersCx, allclustersCy,allclustersCz,2,'filled')
title('Overlay all clusters');
xlabel('x (nm)');
ylabel('y (nm)');
zlabel('z (nm)');
axis equal
% axis([-100 100 -100 100])

binCenters = -200:10:200;
x=transpose(hist(allclustersCx,binCenters)); 
x2=transpose(hist(allclustersCy,binCenters)); 
x3=[x/sum(x)];
x4=[x2/sum(x2)];

subplot(3,2,2)

bar(binCenters, x3)
title('Histogram over x');
xlabel('x (nm)');
ylabel('norm counts');
axis([-100 100 0 0.4])

subplot(3,2,4)
bar(binCenters, x4)
title('Histogram over y');
xlabel('y (nm)');
ylabel('norm counts');
axis([-100 100 0 0.4])

% binCenters = -300:10:300;
x5=transpose(hist(allclustersCz,binCenters));

subplot(3,2,6)
bar(binCenters, x5/sum(x5))
title('Histogram over z');
xlabel('z (nm)');
ylabel('norm counts');
axis([-200 200 0 0.2])


%% Calculate gaussian PDF of x and y dimensions, overlay with histogram
clc
%%%%%%%%%%%%
% Fit the lateral spread of localizations to normal distribution
% Average for all molecules
%%%%%%%%%%%%

figure('Position',[100 500 1000 300],'name','PDF of x and y radius (Overlay)')

% create normal distribution

pdx=fitdist(allclustersCx,'normal')
pdy=fitdist(allclustersCy,'normal')
pdz=fitdist(allclustersCz,'normal')

y  = pdf(pdx,allclustersCx);
y2 = pdf(pdy,allclustersCy);
y3 = pdf(pdz,allclustersCz);

binCenters = -200:10:200;
x  = transpose(hist(allclustersCx,binCenters)); 
x2 = transpose(hist(allclustersCy,binCenters)); 
x5 = transpose(hist(allclustersCz,binCenters)); 
x3 = [x/sum(x)];
x4 = [x2/sum(x2)];
x6 = [x5/sum(x5)];

subplot(1,3,1)
bar(binCenters,x3/max(x3));
hold on
scatter(allclustersCx,y/max(y),1,'red')
axis([-200 200 0 1])
title('PDF x dimension');
xlabel('radius (nm)');
ylabel('pdf');
hold off

clear f x

subplot(1,3,2)
bar(binCenters,x4/max(x4));
hold on
scatter(allclustersCy,y2/max(y2),1,'red')
axis([-200 200 0 1])
title('PDF y dimension');
xlabel('radius (nm)');
ylabel('pdf');
hold on

subplot(1,3,3)
bar(binCenters,x6/max(x6));
hold on
scatter(allclustersCz,y3/max(y3),1,'red')
axis([-200 200 0 1])
title('PDF y dimension');
xlabel('radius (nm)');
ylabel('pdf');
hold on
