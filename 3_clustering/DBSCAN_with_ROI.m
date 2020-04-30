%% DBSCAN clustering with manual ROI

% Author: Christian Sieben, EPFL 
% sieben.christian@gmail.com
% April 2020

% Input:    1. localization table
%           2. xCol
%           3. yCol
%           4. pxlsize for the ROI selection
%           5. k, number of locs within Eps
%           6. Eps, radius
%           
% Output:   1. clustered locs
%           2. FiC, fraction in cluster
%           3. CpA, clusters per area
%           4. LpC, locs per cluster
 
function [subset, FiC, CpA, LpC] = DBSCAN_with_ROI(locs, xCol, yCol, pxlsize, k, Eps);

% pxlsize = 200;

heigth  = round((max(locs(:,yCol))-min(locs(:,yCol)))/pxlsize);
width   = round((max(locs(:,xCol))-min(locs(:,xCol)))/pxlsize);
im      = hist3([locs(:,xCol),locs(:,yCol)],[width heigth]); % heigth x width

% Select rectangle

rect = []; 

figure('Position',[100 200 600 600])
f = imagesc(imrotate(im,90),[0 20]);axis square;
colormap('hot'); colorbar;

rect = getrect; % [xmin ymin width height]

fprintf('\n -- ROI selected --\n')

xmin = min(locs(:,xCol)) + rect(1,1)*pxlsize;
ymin = max(locs(:,yCol)) - rect(1,2)*pxlsize - (rect(1,4)*pxlsize) ;
xmax = xmin + rect(1,3)*pxlsize;
ymax = ymin + rect(1,4)*pxlsize;

target      = find(locs(:,xCol)>xmin & locs(:,xCol)<xmax & locs(:,yCol)>ymin & locs(:,yCol)<ymax);
locs_ROI    = locs(target,1:end);

% figure
% scatter(locs_ROI(:,xCol),locs_ROI(:,yCol),'.')

fprintf('\n -- Plotted selected ROI  --\n')


%% Cluster DBSCAN

% k   = 5;                                                    % minimum number of neighbors within Eps
% Eps = 50;                                                   % minimum distance between points, nm

fprintf('\n -- Parameters selected --\n')

tic
[class,type]    = DBSCAN(locs_ROI(:,xCol:yCol),k,Eps);      % uses parameters specified at input
class2          = transpose(class);                         % class - vector specifying assignment of the i-th object to certain cluster (m,1)
type2           = transpose(type);                          % (core: 1, border: 0, outlier: -1)

fprintf(' -- DBSCAN computed in %f sec -- \n',toc)

%% Find core points --> type = 1 or 0 and define them as subset

coreBorder = find(type2 >= 0);

subset = [];
subset = locs_ROI(coreBorder,1:end);
subset(:,end+1) = class2(coreBorder);


figure('Position',[100 500 600 300])
subplot(1,2,1)
scatter(locs_ROI(:,xCol),locs_ROI(:,yCol),1);
title('input points');xlabel('x [nm]');ylabel('y [nm]');
axis square
axis([min(locs_ROI(:,xCol)) max(locs_ROI(:,xCol)) min(locs_ROI(:,yCol)) max(locs_ROI(:,yCol))])
box on


subplot(1,2,2)
scatter(subset(:,xCol),subset(:,yCol),1,mod(subset(:,end),10))
title('identified clusters');xlabel('x [nm]');ylabel('y [nm]');
axis square
axis([min(locs_ROI(:,xCol)) max(locs_ROI(:,xCol)) min(locs_ROI(:,yCol)) max(locs_ROI(:,yCol))])
box on

clusterlength = [];
for i = 1:max(subset(:,end)); 
    target          = find(subset(:,end)==i);
    clusterlength(i,1)   = length(target);    
end


Area  = (max(locs_ROI(:,xCol))-min(locs_ROI(:,xCol)))*(max(locs_ROI(:,yCol))-min(locs_ROI(:,yCol)));
FiC   = length(subset)/length(locs);  % fraction in cluster  
CpA   = max(class)/(Area/1e6);        % clusters per area/ µm2
LpC   = mean(clusterlength);          % locs per cluster

end

