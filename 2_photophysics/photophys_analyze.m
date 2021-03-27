%% Analyze dye photophysics from SMLM data (processing part)
% 
% Author: Christian Sieben, EPFL 
% sieben.christian@gmail.com
% April 2020

% 1.    Load Tracks
% 2.    Calculate Track length
% 3.1   Calculate Track Spread (i.e. sigma)
% 3.2   Calculate gaussian PDF of x and y dimensions (overlay all molecules)
% 3.3   Calculate gaussian PDF of x and y dimensions (individual molecules)
% 4.    Find Gaps and measure length --> Dark time

clear, clc, close all
SMLM_tutorial_main = '/Users/christian/Documents/Arbeit/MatLab/SMLM_tutorial';

%% 1. Load Tracks

% Output from photophys_head

% To get the merging parameters, use output: Tracks_GapMax_
% The array should have the form:

% res1 = x
% res2 = y
% res3 = photons
% res4 = uncertainty
% res5 = frame
% res6 = track ID

load('A647_COT_1200mW_10ms_3_MMStack_1_Localizations_DC_Z_tracks_GapMax.mat')
all_tracks = res;

%% 2. Calculate Track length

%%%%%%%%%%%%
% On-time distribution (On-time)
%%%%%%%%%%%%
% On time, not the number of blinks
% This is the time (frames) the molecule is visible
% To get the number of blinks, the data need to be merged with 0 gap
%%%%%%%%%%%%


tracklength=[];

for index = 1:max(all_tracks(:,6))
    
    track = find(all_tracks(:,6)==index);
    tracklength = cat(1,tracklength,length(track));
    
    clear track
end

% Plot Track length histogram

binCenters = 0:5:50;
x = transpose(hist(tracklength,binCenters)); 

figure('Position',[100 500 300 300],'name','Track length histogram') 
bar(binCenters,x/sum(x));hold on;
axis([-5 50 0 0.9]);
title(['Mean = ', num2str(mean(tracklength))]);
% xlabel('on time (frames)');
xlabel('number of blinks');
ylabel('norm counts'); 
axis square;


%% 3.1 Calculate Track Spread (i.e. sigma)
%  3.1 can be used to make a figure visualizing sigma

%%%%%%%%%%%%
% The lateral spread of localizations (i.e. sigma)
%%%%%%%%%%%%

clusterx=[];
clustery=[];

allclustersCx=[];
allclustersCy=[];
                                               
for index=1:max(all_tracks(:,6)); % for all tracks
    
    vx=find(all_tracks(:,6)==index);
    
    clusterx=[];
    clustery=[];
    
        if  length(vx)>1;                                                             % if nan, copy  frame number in new

            clusterx=all_tracks(vx,1);
            clustery=all_tracks(vx,2);
            
            clusterxC=sum(clusterx)/length(clusterx);
            clusteryC=sum(clustery)/length(clustery);
            clusterx=clusterx-clusterxC;
            clustery=clustery-clusteryC;
            
            if clusterx>50; 
            else
            

            allclustersCx=vertcat(allclustersCx,clusterx);
            allclustersCy=vertcat(allclustersCy,clustery);
            
            end
        else end

end

clear max

figure('Position',[100 100 500 500],'name','XY spread of all molecules')

subplot(2,2,1)
scatter(allclustersCx, allclustersCy,'k.')
title('overlay all points');
box on; axis square; 
xlabel('x (nm)');
ylabel('y (nm)');
axis([-100 100 -100 100])

subplot(2,2,2)

hist3([allclustersCx, allclustersCy],[50 50])
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
title('All clusters');
xlabel('x (nm)');
ylabel('y (nm)');
colormap('hot')

subplot(2,2,3)
histogram(allclustersCx,40,'Normalization','probability');
title('Histogram over x');
xlabel('x (nm)');
ylabel('norm counts');
axis([-100 100 0 0.3])

subplot(2,2,4)
histogram(allclustersCy,40,'Normalization','probability');
title('Histogram over y');
xlabel('y (nm)');
ylabel('norm counts');
axis([-100 100 0 0.3])

% make another figure

figure('Position',[600 100 500 500]);

ah1 = subplot(4, 4, [3 4 7 8]);
scatter(allclustersCx, allclustersCy,1,'black');
set(gca,'YTicklabel',[],'XTicklabel',[]);
% xlabel('x (nm)','FontSize',10);
% ylabel('y (nm)','FontSize',10);
axis([-100 100 -100 100])
box on 
axis square

subplot(4, 4, [2 6]);
histogram(allclustersCy,40,'Normalization','probability');
xlabel('y (nm)','FontSize',10);
ylabel('norm counts','FontSize',10);
axis([-100 100 0 0.3])
view(-90, -90)

subplot(4, 4, [11 12]);
histogram(allclustersCx,40,'Normalization','probability');
xlabel('x (nm)','FontSize',10);
ylabel('norm counts','FontSize',10);
axis([-100 100 0 0.3])
view(180, 90)

%% 3.2 Calculate gaussian PDF of x and y dimensions (overlay all molecules)
%  3.2 determine average sigma all molecules 

%%%%%%%%%%%%
% Fit the lateral spread of localizations to normal distribution
% Average for all molecules
%%%%%%%%%%%%

figure('Position',[100 500 700 300],'name','PDF of x and y radius (Overlay)')

% create normal distribution

pdx = fitdist(allclustersCx,'normal')
pdy = fitdist(allclustersCy,'normal')

y   = pdf(pdx,allclustersCx);
y2  = pdf(pdy,allclustersCy);


binCenters = -200:10:200;
x=transpose(hist(allclustersCx,binCenters)); 
x2=transpose(hist(allclustersCy,binCenters)); 
x3=[x/sum(x)];
x4=[x2/sum(x2)];

subplot(1,2,1)
bar(binCenters,x3/max(x3));
hold on
scatter(allclustersCx,y/max(y),1,'red')
axis([-100 100 0 1.1])
box on; axis square;
title('PDF x dimension');
xlabel('x (nm)','FontSize',12);
ylabel('pdf','FontSize',12);
hold off
legend(['\sigma _x = ' num2str(pdx.sigma)], 'FontSize',12);

subplot(1,2,2)
bar(binCenters,x4/max(x4));
hold on
scatter(allclustersCy,y2/max(y2),1,'red')
axis([-100 100 0 1.1])
box on; axis square;
title('PDF y dimension');
xlabel('x (nm)','FontSize',12);
ylabel('pdf','FontSize',12);
hold on
legend(['\sigma _y = ' num2str(pdy.sigma)], 'FontSize',12);

%% 3.3 Calculate gaussian PDF of x and y dimensions (individual molecules)
%  3.3 determine average sigma of individial molecules 

%%%%%%%%%%%%
% Fit the lateral spread of localizations to normal distribution
% Per individual molecule
%%%%%%%%%%%%

clusterx=[];
clustery=[];

loc_prec_x=[];
loc_prec_y=[];
                                               
for index=1:max(all_tracks(:,6)); 
    
    vx=find(all_tracks(:,6)==index);
    
    clusterx=[];
    clustery=[];
    
        if  length(vx)>10;                                                             % if nan, copy  frame number in new

            clusterx = all_tracks(vx,1);
            clustery = all_tracks(vx,2);
            
            clusterxC=sum(clusterx)/length(clusterx);
            clusteryC=sum(clustery)/length(clustery);
            clusterx=clusterx-clusterxC;
            clustery=clustery-clusteryC;
            
            pdx=fitdist(clusterx,'normal');
            pdy=fitdist(clustery,'normal');
            
            loc_prec_x = vertcat(loc_prec_x, pdx.sigma);
            loc_prec_y = vertcat(loc_prec_y, pdy.sigma);


        else end

end
clc
Mean_loc_prec_x = fitdist(loc_prec_x,'normal')
Mean_loc_prec_y = fitdist(loc_prec_y,'normal')

%% 4. Find Gaps and measure length --> Dark time

%%%%%%%%%%%%
% Dark time distribution (Off-time)
%%%%%%%%%%%%

allgaps = [];

for i = 1:max(all_tracks(:,6)); 
    
    vx = find(all_tracks(:,6)==i);
    
    track = all_tracks(vx,5)-min(all_tracks(vx,5))+1;
    
    gaps=[];
    gaps(1,1)=1;
   
    for j=2:length(track);
        
    gaps(j,1)=track(j)-track(j-1);
        
    end
    
    gaps=gaps-1;
    gaps=nonzeros(gaps)+1;
  
    allgaps=vertcat(allgaps,gaps);
    
end

binCenters = 0:100:10000;
x = transpose(hist(allgaps,binCenters)); 


figure('Position',[100 500 700 300],'name','Dark time distribution')
subplot(1,2,1)
bar(binCenters,x/max(x));
ax = gca;
ax.YScale = 'log'
title(['Mean = ', num2str(mean(allgaps))]);
xlabel('Off time [frames]');
ylabel('norm counts');
hold on


[CdfY,CdfX] = ecdf(allgaps,'Function','cdf'); 
subplot(1,2,2)
plot(CdfX,CdfY,'LineWidth',3);hold on;
scatter(CdfX(max(find(CdfY < 0.99))),1,'ro','filled');
title(['99 % at ' num2str(CdfX(max(find(CdfY < 0.99)))) ' frames']);
xlabel('Off time [frames]');
ylabel('cumulative probability');
% save('gap.mat','gapfit','ft','gof','binCenters','x')
    
