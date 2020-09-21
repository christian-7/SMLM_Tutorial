% Function to manually select a single rectangular ROI

% Author: Christian Sieben, EPFL 
% sieben.christian@gmail.com
% January 2020

%% 1. Load Localizations

function [locs_ROI] = manualROI(locs,xCol,yCol,pxlsize);
%% 

heigth  = round((max(locs(:,yCol))-min(locs(:,yCol)))/pxlsize);
width   = round((max(locs(:,xCol))-min(locs(:,xCol)))/pxlsize);

figure('Position',[650 400 500 500],'Name','Select ROI')
im=hist3([locs(:,xCol),locs(:,yCol)],[width heigth]); % heigth x width
imagesc(imrotate(im,90),[0 10]);
colormap('hot');axis square; box on;
xlabel('x pxl'); ylabel('y pxl');

rect = getrect; % rect = [xmin ymin width height];

close all

fprintf('\n -- ROI selected --\n')

% Select ROI 
close all

xmin = min(locs(:,xCol))+ rect(:,1)*pxlsize;
ymin = max(locs(:,yCol)) - rect(:,2)*pxlsize - (rect(:,4)*pxlsize);
xmax = xmin + (rect(:,3)* pxlsize);
ymax = ymin + rect(:,4) * pxlsize;

vx      = find(locs(:,xCol)>xmin & locs(:,xCol)<xmax & locs(:,yCol)>ymin & locs(:,yCol)<ymax);
locs_ROI = locs(vx,1:end);

figure('Position',[1200 400 500 500],'Name','Selected ROI')
scatter(locs_ROI(:,xCol)/1e3,locs_ROI(:,yCol)/1e3,'.')
%axis square; 
box on; xlabel('x \mu m'); ylabel('y \mu m');

fprintf('\n -- Plotted selected ROI  --\n')

end
