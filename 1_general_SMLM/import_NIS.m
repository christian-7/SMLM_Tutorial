%% Import data from NIS text file

function [locs] = import_NIS();

[FileName1, PathName1] = uigetfile({'*.csv', 'Tables (*.csv)';'*.txt','Text (*.txt)'}, 'Choose the locs to open.'); % Open Image File

cd(PathName1);

% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 26);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["ChannelName", "X", "Y", "Xc", "Yc", "Height", "Area", "Width", "Phi", "Ax", "BG", "I", "Frame", "Length", "Link", "Valid", "Z", "Zc", "Photons", "LateralLocalizationAccuracy", "Xw", "Yw", "Xwc", "Ywc", "Zw", "Zwc"];
opts.VariableTypes = ["categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
locs = readtable(FileName1, opts);

clear opts

%% 

clear, clc

[FileName1, PathName1] = uigetfile({'*.csv', 'Tables (*.csv)';'*.txt','Text (*.txt)'}, 'Choose the locs to open.'); % Open Image File

cd(PathName1);

clc
opts = delimitedTextImportOptions("NumVariables", 26);
opts.Delimiter = "\t";
opts.DataLines = [2, Inf];
opts.ReadVariableNames = "true"

locs = readtable(FileName1, opts);


