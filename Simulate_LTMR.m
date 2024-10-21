% Simulate GLMs (once fitted)
clear variables; clc
plotFlag = 1; saveFlag = 1; fid = pwd; 
addpath('Functions')
addpath('Filters')
addpath('Stimulus')

% Settings 
Neuron = 'Neuron2';   % Choose Neuron 1 to 4 = Cluster 1 to 4 Neuron
runs = 5;             % Number of trials to simulate
ind_plot = 0;         % Plot Figure (0), No Plot (1)

% Load Filters:
path1 = append('LTMR_glmfit-',Neuron,'.mat');
load(path1);
% Load Stimulus (10Hz - 50Hz - 100Hz at 4.5V):
Test1 = load('10Hz_4.5V.mat'); 
Test2 = load('50Hz_4.5V.mat'); 
Test3 = load('100Hz_4.5V.mat'); 
Signal = [Test1.sig;Test2.sig;Test3.sig]; 

% Table of Basis Fxn Parameters:
Kcel = struct2cell(kbasprs); si = size(Kcel); Kcel = reshape(Kcel,si(2),si(1));
Kmat = cell2mat(Kcel); BasisTab(1,:) = array2table(Kmat);
Hcel = struct2cell(ihbasprs);si = size(Hcel); Hcel = reshape(Hcel,si(2),si(1));
Hmat = cell2mat(Hcel); BasisTab(2,:) = array2table(Hmat);
BasisTab.Properties.VariableNames = {'Lgth','# Eye','# BF','B','Peak 1','Peak N','Refr'};
BasisTab.Properties.RowNames = {'Stim Filt','Spike Filt'};
BasisTab2 = BasisTab(:,[1,3,4,7]);

% Stimulate GLM:
[y, stimcurr, hcurr, r] = simulate_glm(Signal,dt,k,h,dc,runs,[],BasisTab2,softRect,ind_plot);

