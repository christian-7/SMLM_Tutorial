clear, clc, close all

% Use calibrate3D_GUI to generate a PSF calibration file

%% Initialize and scan input directory

fitting_dist = 'X:\Christian\spline_fitting_distribution';                                  % define path to fitting distribution

cd(fitting_dist);

image_stacks        = 'E:\to_analyze\2020-03-04_CS_MTs_Julie\image_stacks';            % Path to MM tiff stacks
output_folder       = 'E:\to_analyze\2020-03-04_CS_MTs_Julie\locResults';              % Path to locResults

diary([output_folder '\logFile.txt']);

input = dir(image_stacks); input_dir = [];

for i = 3:size(input,1);
    
   input_dir{i-2,1}  = [image_stacks '\' input(i).name];

end
                     
%% Initialize the fitter

path_splineFit     = [fitting_dist '\fit3Dcspline\']; % must be a local folder
addpath(genpath(path_splineFit));

% For Miji on Win, run this before runnig starting miji
cd(path_splineFit);

javaaddpath '\ImageJ\ij.jar';
javaaddpath '\ImageJ\mij.jar';

myMiji(true,'ImageJ');
                   
fprintf('\n -- Initialized -- \n'); 

for FOV = 1:size(input_dir,1); % loop through the different image stacks as indexed above
    
clearvars -except input_dir output_folder FOV  fitting_dist

fprintf(['\n -- Starting FOV' num2str(FOV) '  -- \n']); 
    
%% Define input and output parameters

input_folder = input_dir{FOV};

calib_file         = [fitting_dist '\PSF_calibration\MT_Julie_Ch642_ROI17_Filter2_3Dcorr.mat'];          % needs to be generated using the calibratd3D_GUI, described in the tutorial
variance_map       = [fitting_dist '\camera_calibration\prime_alice\varOffset.mat'];        % variance map of the CMOS camera

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(output_folder);
mkdir('temp'); 
temp_output_folder = [output_folder '\temp']; % temporary output, will be deleted at the end

%% Index the input folder

cd(input_folder);

metadata_file = dir(sprintf('*metadata.txt'));

image_files = dir(sprintf('*ome.tif'));

fprintf('\n -- Indexed input folder -- \n');     

%% For sCMOS cameras, load respective variance map and select ROI

% Open metadata and extract the ROI coordinates

fid = fopen(metadata_file.name, 'rt');
s   = textscan(fid, '%s',100, 'delimiter', '\n'); % load only the first 100 lines
fid = fclose(fid);

fprintf('\n -- Metadata file loaded -- \n');   

row = strfind(s{1}, '"ROI":'); roi = {};

for i = 1:length(s{1,1});
    
    if row{i,1} == 1;
    roi = s{1,1}(i+1,1);
    end
    
end

% temp    = regexp(roi,'\"','split');
% roi     = regexp(temp{1,1}{1,4},',','split');

% MM2 beta

temp1       = regexp(roi,'\"','split');
temp2       = regexp(temp1{1,1}{1,4},',','split');
temp3       = regexp(temp2{1,1},'[','split');
temp4       = regexp(temp2{1,4},']','split');

roi = [];
roi(1,1) = str2num(temp3{1,2});
roi(2,1) = str2num(temp2{1,2});
roi(3,1) = str2num(temp2{1,3});
roi(4,1) = str2num(temp4{1,1});


% Extract ROI from camera variance map
% ROI format (x,y,widht,height) --> XY of upper left corner

var_Alice    = load(variance_map); 

var_map = [];
var_map = var_Alice.varOffset(roi(1, 1):roi(1, 1)+roi(3, 1)-1,... 
                              roi(2, 1):roi(2, 1)+roi(4, 1)-1);

                          
save([input_folder '\var_map.mat'],'var_map');                          
                
fprintf('\n -- Variance map selected -- \n');    

%% Loop the fitting routine through the image folder

% This assumes multi-stack image files per FOV as MM restricts the size of
% each tiff stack

clc; startLoop = tic;

for i = 1:size(image_files,1);
    
image_name  = image_files(i).name;
base        = regexp(image_name,'\.','split');

% input_file  = ['path=[' input_folder '\' image_name ']']; % if using normal stack import
input_file  = ['open=[' input_folder '\' image_name ']'];   % if using virtual stack import
output_file = [temp_output_folder '\' base{1} '_Localizations.csv'];

image_files(i).output = output_file; % Save the output path

% MIJ.run('Open...', input_file); % Open tiff stack

MIJ.run('TIFF Virtual Stack...', input_file); % Open virtual tiff stack

startIteration = tic;

p                   = {};
p.imagefile         = '';           % opened in Miji
p.calfile           = calib_file;   % if empty 2D Gaussian will be used
p.offset            = 163;          % in ADU, for Alice: 163
p.conversion        = 0.48;         % e/ADU,  for Alice: 0.48
p.previewframe      = false;
p.peakfilter        = 2;            % filter size (pixel): 1.2/2
p.peakcutoff        = 10;           % photons (this need to be tested using the simpleFitter_GUI) 10/3
p.roifit            = 17;           % ROI size (pixel)
p.bidirectional     = false;        % true for 2D
p.mirror            = false;
p.status            = '';
p.outputfile        = output_file;
p.outputformat      = 'csv';
p.pixelsize         = 106;      % Anna: 110 / Alice: 106
p.loader            = 3;        % {'1 - simple tif','2 - ome loader','3 - ImageJ'}
p.mij               = MIJ;
p.backgroundmode    = 'Difference of Gaussians (fast)';
p.preview           = false;
p.isscmos           = true;
p.scmosfile         = [input_folder '\var_map.mat'];;
p.mirror            = false;

fprintf('\n -- Starting localization ...  \n'); 

simplefitter_cspline(p);

endIteration = toc(startIteration);

fprintf(['\n -- Finished processing substack ' num2str(i) ' of ' num2str(size(image_files,1)) ' in ' num2str(endIteration/60) ' min -- \n']); 

MIJ.run('Close All');

end

endLoop = toc(startLoop);

fprintf(['\n -- Total processing time = ' num2str(endLoop/60) ' min -- \n']); 

%% Combine files and save as single localization file

cd(temp_output_folder); all_locs = [];frames = [];

% Load Header

file      = fopen(image_files(1).output);
line      = fgetl(file);
h         = regexp( line, ',', 'split' );
frameCol  =  strmatch('frame',h);
fclose('all');

for i = 1:size(image_files,1);

    locs     = [];
    locs     = dlmread(image_files(i).output,',',1,0);
    all_locs = vertcat(all_locs,locs);
    
    if isempty(frames)==1;
        frames   = vertcat(frames,locs(:,frameCol));
    else
        frames   = vertcat(frames,locs(:,frameCol)+max(frames));
    end
    
    clc
    fprintf(['\n -- Processed locfile ' num2str(i) ' of ' num2str(size(image_files,1)) '-- \n']); 
    
end

all_locs(:,frameCol) = frames;

delete *.csv
new_name_temp   = regexp(base{1},'_','split');
new_name = [sprintf('%s_',new_name_temp{1:find(strcmp(new_name_temp,'MMStack'))-1}) 'locs.csv'];
% new_name = [sprintf('%s_',new_name_temp{1:end-find(strcmp(new_name_temp,'MMStack'))}),new_name_temp{end-(find(strcmp(new_name_temp,'MMStack'))-1)} '_locs.csv'];
% new_name        = [new_name_temp{1,1} '_' new_name_temp{1,2} '_' new_name_temp{1,3} '_FOV_' new_name_temp{1,4} '_locs.csv'];

cd(output_folder);
rmdir temp

fileID = fopen(new_name,'w');
fprintf(fileID,[[line] ' \n']);
dlmwrite(new_name,all_locs,'-append');
fclose('all');

end

diary off
cd(output_folder);
clc
fprintf('\n -- Finished Processing -- \n'); 