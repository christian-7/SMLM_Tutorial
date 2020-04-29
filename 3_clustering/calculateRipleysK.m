function [K, K_rand] = calculateRipleysK(locs, xCol, yCol, ROI, maxk, stepk);

wb      = waitbar(0,'Computing Ripley''s K-function...');
% ROI     = [xmin xmax ymin ymax];
% maxk    = 100;    % maximum distance
% stepk   = 10;     % step size
% mink    = 1;      % minimum distance

index   = 1;

for k = stepk:stepk:maxk;

waitbar(k/maxk,wb);
    
K(index,2)  = RipleysK(locs(:,xCol:yCol),k,ROI,0);
K(index,1)  = k;

index       = index+1;

end

close (wb); 

% Compare with random dataset

r1 = min(locs(:,xCol))+(max(locs(:,xCol))-min(locs(:,xCol))).*rand(length(locs),1);
r2 = min(locs(:,yCol))+(max(locs(:,yCol))-min(locs(:,yCol))).*rand(length(locs),1);
locs_rand = [r1,r2];

wb      = waitbar(0,'Computing Ripley''s K-function for random dataset...');

K_rand = []; index   = 1;

for k = stepk:stepk:maxk;

waitbar(k/maxk,wb);
    
K_rand(index,2) = RipleysK(locs_rand,k,ROI,0);
K_rand(index,1) = k;

index=index+1;

end

close (wb); 



end
