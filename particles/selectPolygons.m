%% Select multiple polygons

function [particles] = selectPolygons(locs,xCol,yCol);

selection = []; particles = {};

figure('Position',[100 200 600 600])

f = scatter(locs(:,xCol), locs(:,yCol),'k.'); count = 1;

while isvalid(f)

  try  h = impoly();

       nodes = getPosition(h);
       selected_indices = inpolygon(locs(:,xCol),locs(:,yCol), nodes(:,1), nodes(:,2)); 
       particles{count,1} = locs(selected_indices,1:end);
       count = count +1;
       
  catch continue
  end

end

% See if it works

% figure
% scatter(locs(:,x), locs(:,y),'k.');hold on;
% for i = 1:size(particles,1);
% scatter(particles{i, 1}(:,x),particles{i, 1}(:,y),'ro');hold on;
% end

end