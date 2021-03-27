% 1. simulates a random walk according to the initialized parameters and adds confined areas with D_conf
% 2. calculates the confinement index according to Simson, et al. (1995), Biophysical Journal, 69(September), 989?993. 
% 
% 
% D-diffusion coefficient
% D_conf-diffusion coefficient in confined area
% dt-time step
% 
% Date:     22/06/16
% Author:   Christian Sieben

%% 1. Initialize parameters for random walk and add circles

clear all, close all, clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D           = 0.041;                        % diffusion constant, um2/s 
D_conf      = D/100;                        % diffusion constant in confined areas, ?m2/s
dwell_time  = 100;                          % time the particle stays in confined area

start       = [0 0];                        % starting coordinates
num_steps   = 400;                          % number of steps
                    
dt          = 0.5;                          % time step, seconds
dx          = 1;                            % pixel size, um

segment     = 30;                           % Sm, segment length in frames

%% 2. Generate and plot random walk trajectory
close all
% Calculate random diffusion -->tracks.m

k       = sqrt(2 * D * dt);             % Stokes/Einstein Law, sigma of the Normal Dist
k_conf  = sqrt(2 * D_conf * dt);        % Stokes/Einstein Law, sigma of the Normal Dist
pos     = zeros(num_steps,2);
time    = (0 : num_steps-1)' * dt;

well1   = 50; well2 = 300;             % time after which particle falls in well

for n=2:num_steps;
    
    if n>=well1 & n<=well1+dwell_time | n>=well2 & n<=well2+dwell_time
     
    dX = k_conf * randn(1, 2);          % calculate step from normal distribution   
    pos(n,1) = pos(n-1,1)+ dX(1);  
    pos(n,2) = pos(n-1,2)+ dX(2); 
    
    clear dX
    
    else
    
    dX = k * randn(1, 2);    
    pos(n,1) = pos(n-1,1)+ dX(1);  
    pos(n,2) = pos(n-1,2)+ dX(2); 
    
    clear dX
    
    end
    
end

pos(:,3) = time;

% Plot the simulated trajectory

figure('Position',[100 100 500 500], 'name','Simulated Trajectory')
line(pos(:,1)*dx,pos(:,2)*dx,'Color',[0 0 0]);hold on;
% % line(pos(well1:well1+dwell_time,1)*dx,pos(well1:well1+dwell_time,2)*dx,'Color',[1 0 0]);hold on;
% % line(pos(well2:well2+dwell_time,1)*dx,pos(well2:well2+dwell_time,2)*dx,'Color',[1 0 0]);hold on;
axis([min(pos(:,1))-1 max(pos(:,1))+1 min(pos(:,2))-1 max(pos(:,2))+1])
circle1 = viscircles([pos(well1+dwell_time/2,1) pos(well1+dwell_time/2,2)], 0.5);
circle2 = viscircles([pos(well2+dwell_time/2,1) pos(well2+dwell_time/2,2)], 0.5);

title('XY plot trajectory');
xlabel('x (\mum)','FontSize',12);
ylabel('y (\mum)','FontSize',12);
box on; axis square; axis equal

%% Optional: Generate an animated plot

x = pos(:,1);
y = pos(:,2);
fid = figure;
hold on
axis([min(x) max(x) min(y) max(y)]); 
title('XY plot trajectory');
xlabel('x (\mum)','FontSize',12);
ylabel('y (\mum)','FontSize',12);
box on
axis square

circle1 = viscircles([pos(well1+dwell_time/2,1) pos(well1+dwell_time/2,2)], 0.5);
circle2 = viscircles([pos(well2+dwell_time/2,1) pos(well2+dwell_time/2,2)], 0.5);

% Set up the movie.
writerObj = VideoWriter('out2.avi'); % Name it.
writerObj.FrameRate = 60; % How many frames per second.
open(writerObj); 
 
for i=1:size(y)-1      
    % We just use pause but pretend you have some really complicated thing here...
    pause(0.01);
    figure(fid); % Makes sure you use your desired frame.
    plot(x(i:i+1),y(i:i+1),'Color','k','LineWidth', 1.5);
 
    %if mod(i,4)==0, % Uncomment to take 1 out of every 4 frames.
        frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
        writeVideo(writerObj, frame);
    %end
 
end
hold off
close(writerObj); % Saves the movie.



%% 3. Calculate confinement

% generate variable frame

frame = pos(:,3);             % time step in seconds
frame = frame/dt;             % time step in frames
frame = frame-min(frame)+1;     % starting from 0
% frame(1,1)=1;               % starting from 1

% i = frame --> row
% j = gap; --> column

prob=[];    %zeros(max(frame), 5);
prob2=[];   %zeros(5, 2);
d=[];
vx=[];
vy=[];

c=1;

for i=1:max(frame)-segment;                  % for all frames
    vx=find(frame == i);                     % find frame i
    
    if isempty(vx)==1;                       % if frame does not exist, skip   
    else
        
     
    for j=4:segment;                                    % segment length
          
        vy=find(frame <= (i+j) & frame >= i );          % select segment
        subset(:,1)=pos(vy);                            % define segment as subset
        subset(:,2)=pos(vy,2);
        
        if length(vy)==1;   % if subset is only 1 frame --> distance is 0
                     R=0;
        else    
        
            for  k=2:length(subset);
                 d(k,1)=sqrt(((subset(k,1)-subset(1,1))^2)+((subset(k,2)-subset(1,2))^2));    % calculate the distance to each point in subset from point i  
            end
        R=max(d);                                                      % maximum distance within subset
        prob(i:(i+j),c)=0.2048-2.5117*((D*j)./(R^2));                  % probability within subset
%         prob(c,i:(i+j))=((D*j)./(R^2));

%       prob(c,i:(i+j))=horzcat(prob(c,i:(i+j)),((D*j)./(R^2)));
        c=c+1;  
        clear subset
        end
    
%     c=c+1;    
        
    end
    clear vx vy R d;
    
    end
   
end
clear subset

for l=1:max(frame)
    
prob2(l,1)=l;                           % frame
prob2(l,2)=mean(nonzeros(prob(l,:)));   % this is psi

end


L=[];

for i=1:length(prob2)
    
    if 10.^(prob2(i,2))>0.1
       L(i,1)=0;
       
    else
        
        L(i,1)=((prob2(i,2))*(-1)-1);
        
    end
L(i,2)=i;
    
end

% Plot the confinement index vs time

figure('Position',[400 100 600 300],'name','Result Confinement Calculation')
subplot(2,1,1)
plot(L(:,2)*dt, L(:,1),'Color',[0 0 0]); hold on;
xlabel('time [s]','FontSize',14);
ylabel('probability level','FontSize',14);
axis([0 num_steps*dt 0 max(L(:,1))*1.2])
box on;

subplot(2,1,2)
plot(pos(:,3), pos(:,1)-min(pos(:,1)),'Color',[1 0 0]); hold on;
plot(pos(:,3), pos(:,2)-min(pos(:,2)),'Color',[0 1 0]); hold on;
xlabel('time [s]','FontSize',14);
ylabel('XY displacement','FontSize',14);
axis([0 num_steps*dt 0 max([pos(:,1)-min(pos(:,1));pos(:,2)-min(pos(:,2))])*1.2])
box on;

figure('Position',[100 100 300 300], 'name','Simulated Trajectory')
line(pos(:,1)*dx,pos(:,2)*dx,'Color',[0 0 0]);hold on;
axis([min(pos(:,1))-1 max(pos(:,1))+1 min(pos(:,2))-1 max(pos(:,2))+1])
title('XY scatter trajectory');
xlabel('x (\mum)','FontSize',12);
ylabel('y (\mum)','FontSize',12);
box on;



%% 4. Figures to save

close all

figure('Position',[600 800 600 300], 'name','XY scatter and confinement index L')
h=gcf;
set(h,'PaperOrientation','landscape');

subplot(1,2,1)
line(pos(:,1)*dx,pos(:,2)*dx,'Color',[0 0 0]);hold on;
scatter(pos(:,1)*dx,pos(:,2)*dx,10,pos(:,3),'filled');hold on;
plot(pos(1,1)*dx,pos(1,2)*dx,'*b','MarkerSize',12);hold on;
text(pos(1,1)*dx,pos(1,2)*dx, 'Start');
plot(pos(length(pos),1)*dx,pos(length(pos),2)*dx,'+b','MarkerSize',12);hold on;
text(pos(length(pos),1)*dx,pos(length(pos),2)*dx,'End');hold on;
circle1 = viscircles([pos(well1+dwell_time/2,1) pos(well1+dwell_time/2,2)], 0.5);
circle2 = viscircles([pos(well2+dwell_time/2,1) pos(well2+dwell_time/2,2)], 0.5);
title('XY scatter trajectory');
xlabel('x (\mum)','FontSize',12);
ylabel('y (\mum)','FontSize',12);
box on; axis equal;

subplot(1,2,2)
plot(L(:,2)*dt, L(:,1)); hold on;
scatter(L(:,2)*dt, L(:,1),15,pos(:,3),'filled')     
xlabel('time (s)','FontSize',12);
ylabel('confinement index L','FontSize',12);
title('Confinement index L');





figure('Position',[100 800 500 600], 'name','Compare time with displacement and confinement')
h=gcf;
set(h,'PaperOrientation','portrait');

subplot(4,1,1)
line(frame*dt,pos(:,1));hold on;
scatter(frame*dt,pos(:,1),5,pos(:,3)*dt);hold on;
xlabel('time (s)','FontSize',12);
ylabel('x (\mu m)','FontSize',12);
box on;



subplot(4,1,2)
line(frame*dt,pos(:,2));hold on;
scatter(frame*dt,pos(:,2),5,pos(:,3));hold on;
xlabel('time (s)','FontSize',12);
ylabel('y (\mu m)','FontSize',12);
box on;

subplot(4,1,3)
plot(L(:,2)*dt, L(:,1)); hold on;
scatter(L(:,2)*dt, L(:,1),15,pos(:,3),'filled')     
xlabel('time (s)','FontSize',12);
ylabel('confinement index L','FontSize',12);
box on;

dcum = zeros(num_steps,1);
dcum(1,1)=0;

for k=2:num_steps;
    
    dist(k,1)=sqrt(((pos(k,1)-pos(1,1))^2)+((pos(k,2)-pos(1,2))^2));

    dcum(k)=dist(k)+dcum(k-1);
end


subplot(4,1,4)
plot(frame*dt,dist); hold on;
scatter(frame*dt,dist,2,pos(:,3)*dt);
xlabel('time (s)','FontSize',12);
ylabel('distance from origin (\mu m)','FontSize',12);
box on;

%% 5. Identify wells, measure dwell time and size

threshold = 3; % if confinement index raises above its counted as confined

dwell = [];
dwell_find = find(L(:,1) > threshold); % Find the wells

dwell(:,1) = dwell_find;
dwell(:,2) = pos(dwell_find,1);
dwell(:,3) = pos(dwell_find,2);
dwell(:,1) = dwell(:,1)-min(dwell(:,1))+1;

wells=[]; well_scatter = {};
dwell_time=[];
p = 1;

for o=1:length(dwell)-1;
    
    if dwell(o+1,1)==dwell(o,1)+1;
    
    wells(o,p) = dwell(o,1);
    
    well_scatter{p,1}(o,1) = dwell(o,2);
    well_scatter{p,1}(o,2) = dwell(o,3);

    else
        
    p = p+1;
    wells(o,p)=dwell(o,1);

    end
    
end

[m,n]=size(wells);

for q=1:n;
    
dwell_time(q,1)=length(nonzeros(wells(:,q)))

end

% figure('Position',[500 800 400 400], 'name','Dwell time histogram')
% h=gcf;
% set(h,'PaperOrientation','portrait');
% hist(dwell_time) 

figure('Position',[500 800 400 400], 'name','Dwell time histogram')

for i = 1:length(well_scatter);
    
    x = nonzeros(well_scatter{i,1}(:,1)); x = x - sum(x)/length(x);
    y = nonzeros(well_scatter{i,1}(:,2)); y = y - sum(y)/length(y);
    
    k = convhull(x,y);

    plot(x(k),y(k),'r-');hold on;
    
    
    
end
