%% calculate MSD and D

% Author: Christian Sieben, EPFL 
% sieben.christian@gmail.com
% January 2020

% This script uses the @msdanalyzer class developed by Jean-Yves Tinevez
% download at: https://tinevez.github.io/msdanalyzer/

% Content

% 1. Initialize @msdanalyzer 
% 2. Calculate MSD curves
% 3. Estimating diffusion coefficient

close all
%% 1. Initialize @msdanalyzer 

% initiate units

SPACE_UNITS     = 'um';
TIME_UNITS      = 's';

% initiate msdanalyzer

ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
ma = ma.addAll(tracks);

%% 2. Calculate MSD curves

% Compute MSD and plot individual MSD curves

ma = ma.computeMSD;
ma.msd

figure('Position',[500 600 600 200],'Name',['MSD Analysis'])
subplot(1,2,1)
ma.plotMSD; box on
title('individual tracks')
subplot(1,2,2)
ma.plotMeanMSD(gca, true);box on
title('weighted average curve')

%% 3. Estimating diffusion coefficient
clc

% Mean of all curves

[fo, gof] = ma.fitMeanMSD;
fprintf('\n');

% Individual curves

ma = ma.fitMSD;

good_enough_fit = ma.lfit.r2fit > 0.8;
Dmean = mean( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;
Dstd  =  std( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;

fprintf('Estimation of the diffusion coefficient from linear fit of the MSD curves:\n')
fprintf('D = %.3g +/- %.3g (mean +/- std, N = %d)\n', ...
Dmean, Dstd, sum(good_enough_fit));
