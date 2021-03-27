clear all, close all, clc

% Generate Clusters without noise --> store in variable clusters

a = 5;       % lower XY limit
b = -5;      % upper XY limit

centersx = (b-a).*rand(10,1) + a; % generate centers
centersy = (b-a).*rand(10,1) + a; % generate centers


count=1;
clusterx=[];
clustery=[];

n=200; % number of points per cluster

for count=1:length(centersx);
   
    clustercx = (randn(n,1)/10)+centersx(count);
    clustercy = (randn(n,1)/10)+centersy(count);
    clusterx=[clusterx; clustercx];
    clustery=[clustery; clustercy];
    clear clustercx cluster cy
    count=count+1;
   
end

clusters(:,1)=clusterx;
clusters(:,2)=clustery;

figure
scatter(clusters(:,1),clusters(:,2),'.');
axis square
box on

%% Add noise

noisex = (b-a).*rand(1000,1) + a;
noisey = (b-a).*rand(1000,1) + a;

clusterx=[clusterx; noisex];
clustery=[clustery; noisey];

clear clusters

clusters(:,1)=clusterx;
clusters(:,2)=clustery;

figure
scatter(clusters(:,1),clusters(:,2),'.');
axis square
box on

%% Set up some function. 
% Sine between -2*pi and 2*pi.
x = (10*-pi:0.1:10*pi)'; % Note the transpose.
y = sin(x);
fid = figure;
hold on
% The final plot.
% plot(x,y, '*');
axis([min(x) max(x) min(y) max(y)]); 
title('XY plot trajectory');
xlabel('x (\mum)','FontSize',12);
ylabel('y (\mum)','FontSize',12);
box on
axis square
%% Set up the movie.
writerObj = VideoWriter('out.avi'); % Name it.
writerObj.FrameRate = 60; % How many frames per second.
open(writerObj); 
 
for i=1:size(y)-1      
    % We just use pause but pretend you have some really complicated thing here...
    pause(0.01);
    figure(fid); % Makes sure you use your desired frame.
    plot(x(i:i+1),y(i:i+1),'Color','k','LineWidth', 1.5);
 
    %if mod(i,4)==0, % Uncomment to take 1 out of every 4 frames.
%         frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
        writeVideo(writerObj, frame);
    %end
 
end
hold off
close(writerObj); % Saves the movie.

