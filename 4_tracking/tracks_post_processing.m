%% Post-process sptPALM tracks

% Author: Christian Sieben, EPFL 
% sieben.christian@gmail.com
% January 2020

% Content:
% 
% 1. Identify Mobile/Immobile fraction
% 2. Histogram of D
% 3. Plot Mean MSD and STD MSD for mobile/immobile tracks
%
%% 1. Identify Mobile and Immobile molecules according to D

close all

mobile      = []; immobile    = []; coords      = []; countM      = 1; countIM     = 1;

figure('Position',[100 500 900 450])
subplot(1,2,2);

for index=1:length(ma.lfit.a);
    
    coords = vertcat(coords,ma.tracks{index, 1}(:,2:3)); 
    
    if ma.lfit.a(index,1)>0.05 & ma.lfit.a(index,1)<1 % D > 0.05 --> mobile
        
         mobile(:,countM)=ma.msd{index,1}(:,2); % one track per column
         countM=countM+1;
         plot(ma.tracks{index, 1}(:,2),ma.tracks{index, 1}(:,3),'black');hold on;
 
    elseif  ma.lfit.a(index,1)>0 & ma.lfit.a(index,1)<0.05
            immobile(:,countIM)=ma.msd{index,1}(:,2); 
            countIM=countIM+1;
            plot(ma.tracks{index, 1}(:,2),ma.tracks{index, 1}(:,3),'red');hold on;
        
    else
        
             
    end
    
end

title(['Tracks: ' num2str(countIM/length(ma.lfit.a),2) ' immobile / ' num2str(countM/length(ma.lfit.a),2) ' mobile'],'FontSize',14);
legend('Mobile','Immobile','FontSize',12); xlabel('x [nm]','FontSize',12); ylabel('y [nm]','FontSize',12); axis square

pxlsize = 20; % nm 

heigth  = round((max(coords(:,2))-min(coords(:,2)))*pxlsize);
width   = round((max(coords(:,1))-min(coords(:,1)))*pxlsize);
im      = hist3([coords(:,1),coords(:,2)],[width heigth]); % heigth x width

subplot(1,2,1)
imagesc(imgaussfilt(imrotate(im,90),1),[0 6]);
title('Rendered Tracks','FontSize',14);
xlabel('x [nm]','FontSize',12); ylabel('y [nm]','FontSize',12);axis square
colormap hot

%% 2. Histogram of diffusion constants

count=1;

for index=1:length(ma.lfit.a);
    
    if ma.lfit.a(index,1)>0 & ma.lfit.a(index,1)<1 % D > 0.05 --> mobile
        
        FORHist(count,1)=ma.lfit.a(index,1);
        count=count+1;
    else
    end
end


figure('Position',[500 800 250*1.5 210*1.5],'Name','Dist of D')
h=gcf;
set(h,'PaperOrientation','landscape');

binCenters = 0:0.075:1;
x=transpose(hist(FORHist,binCenters)); 
x=x/sum(x);


bar(binCenters,x,'k');hold on;
axis([-0.05 1 0 0.8]);
xlabel('D [\mum^2/s] ','FontSize',12);
ylabel('norm counts','FontSize',12);

%% 3. Plot Mean MSD and STD MSD for mobile/immobile tracks

close all

meanMobile  = [];
stdMobile   = [];

meanIMMobile    = [];
stdIMMobile     = [];

maxLag = 6; % longest lag time to consider

for index=1:maxLag % length(mobile)
    
meanMobile(index,1) = mean(mobile(index,not(isnan(mobile(index,:)))));
stdMobile(index,1)  = std(mobile(index,not(isnan(mobile(index,:)))));

end

for index=1:maxLag%length(mobile)

meanIMMobile(index,1)=mean(immobile(index,not(isnan(immobile(index,:)))));
stdIMMobile(index,1)=std(immobile(index,not(isnan(immobile(index,:)))));

end

figure('Position', [100, 800, 250*1.5, 210*1.5]);; hold on;

s1 = scatter(ma.msd{1,1}(1:maxLag,1),meanMobile,'sb');hold on;
plot(ma.msd{1,1}(1:maxLag,1),meanMobile,'b')
errorbar(ma.msd{1,1}(1:maxLag,1),meanMobile,stdMobile,'b'); 

s2=scatter(ma.msd{1,1}(1:maxLag,1),meanIMMobile,'or');hold on;
plot(ma.msd{1,1}(1:maxLag,1),meanIMMobile,'r');
errorbar(ma.msd{1,1}(1:maxLag,1),meanIMMobile,stdIMMobile,'r'); 

legend([s1 s2],{'Mobile, D>0.05','Immobile, D<0.05'},'FontSize',12)

xlabel('lag time [s]','FontSize',12);
ylabel('MSD [\mum^2]','FontSize',12)
box on;

figure('Position', [400, 800, 250*1.5, 210*1.5]); hold on;
s1=scatter(ma.msd{1,1}(1:maxLag,1),meanMobile,'sb');hold on;
shadedErrorBar(ma.msd{1,1}(1:maxLag,1),meanMobile,stdMobile,'-b',1);

s2=scatter(ma.msd{1,1}(1:maxLag,1),meanIMMobile,'or');hold on;
shadedErrorBar(ma.msd{1,1}(1:maxLag,1),meanIMMobile,stdIMMobile,'-r',1);
legend([s1 s2],{'Mobile, D>0.05','Immobile, D<0.05'},'FontSize',12)

xlabel('lag time [s]','FontSize',12);
ylabel('MSD [\mum^2]','FontSize',12)
box on;
