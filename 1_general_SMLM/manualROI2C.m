% Function to manually select a single rectangular ROI

% Author: Christian Sieben, EPFL 
% sieben.christian@gmail.com
% January 2020

%% 1. Load Localizations

function [locs_ROIC1, locs_ROIC2] = manualROI2C(locsC1,locsC2,xCol,yCol,pxlsize);
%% 

heigth  = round((max(locsC1(:,yCol))-min(locsC1(:,yCol)))/pxlsize);
width   = round((max(locsC1(:,xCol))-min(locsC1(:,xCol)))/pxlsize);

figure('Position',[650 400 500 500],'Name','Select ROI')
im=hist3([locsC1(:,xCol),locsC1(:,yCol)],[width heigth]); % heigth x width
imagesc(imrotate(im,90),[0 10]);
colormap('hot');axis square; box on;
xlabel('x pxl'); ylabel('y pxl');

rect = getrect; % rect = [xmin ymin width height];

close all

fprintf('\n -- ROI selected --\n')

% Select ROI 
close all

xmin = min(locsC1(:,xCol))+ rect(:,1)*pxlsize;
ymin = max(locsC1(:,yCol)) - rect(:,2)*pxlsize - (rect(:,4)*pxlsize);
xmax = xmin + (rect(:,3)* pxlsize);
ymax = ymin + rect(:,4) * pxlsize;

vx          = find(locsC1(:,xCol)>xmin & locsC1(:,xCol)<xmax & locsC1(:,yCol)>ymin & locsC1(:,yCol)<ymax);
locs_ROIC1  = locsC1(vx,1:end);

vx          = find(locsC2(:,xCol)>xmin & locsC2(:,xCol)<xmax & locsC2(:,yCol)>ymin & locsC2(:,yCol)<ymax);
locs_ROIC2  = locsC2(vx,1:end);

figure('Position',[1200 400 500 500],'Name','Selected ROI')
scatter(locs_ROIC1(:,xCol)/1e3,locs_ROIC1(:,yCol)/1e3,'.');hold on;
scatter(locs_ROIC2(:,xCol)/1e3,locs_ROIC2(:,yCol)/1e3,'.');hold on;
axis square; box on; xlabel('x \mu m'); ylabel('y \mu m');

fprintf('\n -- Plotted selected ROI  --\n')

end
