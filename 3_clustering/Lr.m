%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Juliette Griffie & Christian Sieben
%Latest update: 05/05/2020
%version alpha
%Aim: local linearised Ripley's K density estimator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialisation

radius  = 50;
ROI     = 1000;

% Generated
nb_points = 1000;
data = (ROI+2*radius)*rand(nb_points, 2)-radius;

figure()
scatter(data(:,1),data(:,2),'.')

% Loaded -> (x,y) coordinates or (x,y,z)
%nb_points=size(data,1);

%% Main
% We need to discuss about borders

D = pdist2(data,data);
L = [];

for i= 1:nb_points
    delta = D(i,:);
    delta(delta>radius)=[];
    L = [L; sqrt((ROI*ROI/(pi*nb_points))*(size(delta,2)-1))];
    %L=[L; sqrt((ROI*ROI/(pi*nb_points))*(size(delta,2)-1))-radius];
end

data=[data, L];
L=[];
data(data(:,1)<0,:)=[];
data(data(:,2)<0,:)=[];
data(data(:,1)>1000,:)=[];
data(data(:,2)>1000,:)=[];

figure()
box on
hold on
scatter(data(:,1),data(:,2),30, data(:,3),'filled')

%% FYI L curves

radius = 5:5:200;
ROI    = 1000;

% Generated
nb_points       = 2000;
data            = (ROI+2*max(radius))*rand(nb_points, 2)-max(radius);
data_cropped    = data;
data_cropped(data_cropped(:,1)<0,:) = [];
data_cropped(data_cropped(:,2)<0,:) = [];
data_cropped(data_cropped(:,1)>1000,:) = [];
data_cropped(data_cropped(:,2)>1000,:)=[];

% Loaded -> (x,y) coordinates or (x,y,z)
%nb_points=size(data,1);

D = pdist2(data_cropped,data);
D_vector = D(:);
D_vector(D_vector==0)=[];
D_vector(D_vector>max(radius))=[];
H = histogram(D_vector, size(radius,2));
count = H.Values;

L_curves = [];
L_curves = [L_curves,0];

for i=1:size(radius,2)
    
    L_curves=[L_curves, sqrt((ROI*ROI/(pi*(size(data_cropped,1)^2)))*sum(count(1:i)))-radius(i)];
    
end
radius=[0,radius];

figure()
box on
hold on
plot(radius,L_curves)

