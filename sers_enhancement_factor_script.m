% Sample script to test function sers_enhancement_factor.m
%% Load SERS data
load('SERS_data_example.mat');

%% Check toolbox dependencies
% dependencies.toolboxDependencyAnalysis({'sers_enhancement_factor.m'})

%% Compute SERS enhancement factor
[EF, paramStruct] = sers_enhancement_factor(waveNumber, ramanSpectrum, SERSspectrum,...
    'band', 1107, 'lambda', 570, 'molWeight', 550, 'NA', 0.25, 'rho', 1.26, ...
    'surfArea', 8, 'showPlot', true);

% EOF
